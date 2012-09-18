/* *****************************************************************************
 *
 *   DWSOLVER - a general, parallel implementation of the Dantzig-Wolfe
 *               Decomposition algorithm.
 *
 *   Copyright 2010 United States Government National Aeronautics and Space
 *   Administration (NASA).  No copyright is claimed in the United States under
 *   Title 17, U.S. Code. All Other Rights Reserved.
 *
 *   This program is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   This program is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   In accordance with the GNU General Public License Version 3, 29 June 2007
 *   (GPL V3) Section 7. Additional Terms, Additional Permissions are added as
 *   exceptions to the terms of the GPL V3 for this program.  These additional
 *   terms should have been received with this program in a file entitled
 *   "ADDITIONAL_LICENSE_TERMS".  If a copy was not provided, you may request
 *   one from the contact author listed below.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 *   Contact: Joseph Rios <Joseph.L.Rios@nasa.gov>
 *
 **************************************************************************** */

#include <glpk.h>        /* For all of the main GLPK stuff. */
#include <stdio.h>       /* printf, etc */
#include <stdlib.h>      /* malloc, etc */
#include <pthread.h>     /* For threading. */

#ifdef USE_INTEL_MKL
#include <mkl_spblas.h>  /* Sparse matrix stuff. */
#include <mkl.h>         /* Other Math Kernel Library stuff. */
#endif

#include "dw_blas.h"
#include "dw.h"
#include "dw_subprob.h"
#include "dw_support.h"
#include "dw_phases.h"

/* This is the subproblem thread.  It does all the work of the subproblems.
 * When launched, it sets up some data items (housekeeping stuff mostly).  The
 * subproblem then does an initial solve of the problem it has been assigned.
 * Then it checks if the master is ready (the master may be busy launching other
 * subproblems still or setting up its own data).  Next, it sets up a column
 * translation data structure to aid in passing data back to the master.  At
 * that point, we are ready to enter the iterative loop of communication and
 * solving with the master problem.  There are more comments throughout the
 * code below.s
 */
void* subproblem_thread(void* arg) {
	int i, my_col, master_rows, master_cols, j, ret;
	int num_clients, verbosity;
	int temp, len, max_iterations;
	int id;
	double* my_solution;
	double* y;
	double* d_vector;

#ifdef USE_INTEL_MKL
	static char no_trans  = 'N';          // Transpose flags for matrix/vector ops.
	static char yes_trans = 'T';          //
	static char* matdescra = "GXXC";      // Description of the matrix used in an op.
	static char* matdescra_f = "GXXF";
	static int increment = 1;             // A param for stepping thru vectors.
#endif

	static double double_neg_one = -1.0;

	char* local_buffer = malloc(sizeof(char)*100);

	double* col_zero_vector;

	subprob_struct* my_data = (subprob_struct*) arg;
	glp_prob* lp;

	glp_smcp *simplex_control_params = (glp_smcp*) malloc(sizeof(glp_smcp));

	/* Make some local copies of of things that won't change globally. */
	id = my_data->my_id;
	num_clients = my_data->globals->num_clients;
	verbosity = my_data->globals->verbosity;

	my_data->simplex_control_params = simplex_control_params;

	/* Set up parameters appropriately. */
	glp_init_smcp(simplex_control_params);
	simplex_control_params->presolve = GLP_OFF;
	simplex_control_params->msg_lev  = GLP_MSG_ERR;

	/* Read in assigned problem and solve to get initial feasible solution. */
	pthread_mutex_lock(&glpk_mutex);
	lp = lpx_read_cpxlp(my_data->infile_name);
	pthread_mutex_unlock(&glpk_mutex);
	glp_iocp* int_parm = malloc(sizeof(glp_iocp));
	glp_init_iocp(int_parm);
	int_parm->msg_lev = GLP_MSG_ERR;
	for( i = 1; i <= glp_get_num_cols(lp); i++ ) {
		if( glp_get_col_type(lp, i) == GLP_LO ) {
			glp_set_col_bnds(lp, i, GLP_DB, glp_get_col_lb(lp, i), 3000.0);
		}
	}
	ret = glp_simplex(lp, simplex_control_params);
	my_data->obj= glp_get_obj_val(lp);
	if( my_data->globals->enforce_sub_integrality ) {
		//printf("Problem %d going to integer optimize.\n", my_data->my_id);
		glp_intopt(lp, int_parm);
		my_data->obj = glp_mip_obj_val(lp);
		//printf("Problem %d done integer optimizing.\n", my_data->my_id);
	}
	glp_create_index(lp);

	//dw_printf(IMPORTANCE_DIAG,
	//		"Thread %d: Finished initial solve.  Simplex returns %d.  Checking on master...\n", id, ret);
	my_data->lp = lp;

	/* Check that the master problem is ready for reading.  If not, wait. */
	pthread_mutex_lock(&master_lp_ready_mutex);

	//dw_printf(IMPORTANCE_DIAG, "Thread %d: Do I need to wait for a signal?\n",
	//		id);
	if (!signals->master_lp_ready) {

		//dw_printf(IMPORTANCE_DIAG, "Thread %d: Yep.  Waiting...\n", id);
		pthread_cond_wait(&master_lp_ready_cv, &master_lp_ready_mutex);
		//dw_printf(IMPORTANCE_DIAG,
		//	"Thread %d: Master_lp_ready condition signal received.\n", id);
	}
	else {
		//dw_printf(IMPORTANCE_DIAG, "Thread %d: Nope. Master is ready...\n", id);
	}
	pthread_mutex_unlock(&master_lp_ready_mutex);

	/* Set up data structures that rely on master info. */
	my_col                    = glp_get_num_cols(lp)+1;
	my_data->num_cols         = glp_get_num_cols(lp);
	my_data->num_cols_plus    = my_col;
	master_rows               = my_data->md->rows;
	master_cols               = my_data->md->cols;
	my_solution               = calloc(my_col, sizeof(double));
	y                         = calloc(master_cols, sizeof(double));
	my_data->col_translate    = malloc(sizeof(int)*my_col);
	my_data->c                = calloc(my_col, sizeof(double));
	col_zero_vector           = calloc(my_col, sizeof(double));
	my_data->condensed_x      = malloc(sizeof(double)*my_col);
	d_vector                  = malloc(sizeof(double)*my_col);
	my_data->val              = malloc(sizeof(double)*(D->rows_plus+num_clients));
	my_data->ind              = malloc(sizeof(int)*(D->rows_plus+num_clients));
	my_data->double_vector    = malloc(sizeof(double)*(my_col+num_clients));
	my_data->D_row_coords     = malloc(sizeof(int));
	my_data->D_col_coords     = malloc(sizeof(int));
	my_data->D_nnz		      = 0;
	my_data->D_vals		      = malloc(sizeof(double));
	my_data->current_solution = my_solution;
	my_data->r                = DW_INFINITY;
	my_data->unbounded        = glp_get_status(lp);

	/* The col_translate array is a LUT for translating between my column
	 * indices and the col indices in the master problem.  Specifically, given
	 * column (variable) i in my problem, that column has an index of
	 * col_translate[i] in the master problem.  While we're in this loop, store
	 * the values for my portion of the master 'c' vector.
	 */
	int*    temp_ind = malloc(sizeof(int)*(master_rows+1));
	double* temp_val = malloc(sizeof(double)*(master_rows+1));

	/* Locking the master_mutex is excessive here. The goal is to gather data
	 * related to the master while being assured that data won't change.  Need
	 * to implement a better mechanism here. This implementation serializes this
	 * section of code over all subproblems.
	 */
	pthread_mutex_lock(&master_mutex);
	for( i = 1; i < my_col; i++ ) {
		temp = glp_find_col(original_master_lp, glp_get_col_name(lp, i));
		if( my_data->globals->verbosity >= OUTPUT_ALL && id >= 0 ) {
			printf("  %d:  %s (index = %d) in my problem maps to column %d in original.  col_translate[%d] = %d\n.",
					id, glp_get_col_name(lp, i), i, glp_find_col(original_master_lp, glp_get_col_name(lp, i)), i, temp);
//			dw_printf(IMPORTANCE_DIAG,"    %d: %s ", i, glp_get_col_name(lp, i));
//			dw_printf(IMPORTANCE_DIAG,"in my problem maps to column ");
//			dw_printf(IMPORTANCE_DIAG,"%d ", glp_find_col(original_master_lp,
//							glp_get_col_name(lp, i)));
//			dw_printf(IMPORTANCE_DIAG,"in the original problem.\n");
//			dw_printf(IMPORTANCE_DIAG,"        col_translate[%d] = %d\n", i, temp);
		}

		if( temp == 0 ) {
			dw_printf(IMPORTANCE_VITAL,"YOW, %d %s (%s)\n", i,
					glp_get_col_name(lp, i), my_data->infile_name);
			dw_printf(IMPORTANCE_VITAL,"PROBLEM: Column %d (%s) from file %s ", i,
								glp_get_col_name(lp, i), my_data->infile_name);
			dw_printf(IMPORTANCE_VITAL, "NOT found in the master problem.\n");
			dw_printf(IMPORTANCE_VITAL, "Edit your input files to fix this.\n");
			dw_printf(IMPORTANCE_VITAL,"There will likely be a crash now.\n");
		}

		my_data->col_translate[i] = temp;
		my_data->c[i] = glp_get_obj_coef(original_master_lp, temp);

		len = glp_get_mat_col(original_master_lp, temp, temp_ind, temp_val);
		if( len > 0 ) {
			my_data->D_col_coords = realloc(my_data->D_col_coords, sizeof(int)*(my_data->D_nnz+len));
			if( my_data->D_col_coords == NULL ) {
				dw_printf(IMPORTANCE_VITAL,"realloc() returned null?\n");
				dw_printf(IMPORTANCE_VITAL,"Probably should exit?");
			}
			my_data->D_row_coords = realloc(my_data->D_row_coords, sizeof(int)*(my_data->D_nnz+len));
			if( my_data->D_row_coords == NULL ) {
				dw_printf(IMPORTANCE_VITAL,"realloc() returned null?\n");
				dw_printf(IMPORTANCE_VITAL,"Probably should exit?");
			}
			my_data->D_vals = realloc(my_data->D_vals, sizeof(double)*(my_data->D_nnz+len));
			if( my_data->D_vals == NULL ) {
				dw_printf(IMPORTANCE_VITAL,"realloc() returned null?\n");
				dw_printf(IMPORTANCE_VITAL,"Probably should exit?");
			}
			for( j = 1; j <= len; j++  ) {
				my_data->D_row_coords[my_data->D_nnz + j - 1] = temp_ind[j]-1;
				my_data->D_col_coords[my_data->D_nnz + j - 1] = i-1;
				my_data->D_vals[my_data->D_nnz + j - 1] = temp_val[j];
			}
			my_data->D_nnz += len;
		}
	}
	free(temp_ind);
	free(temp_val);

	pthread_mutex_unlock(&master_mutex);

	/* Create my portion of the original objective coeff's. */
	dw_printf(IMPORTANCE_DIAG,"Thread %d: There are %d cols.\n", id, my_col-1);
	organize_solution(my_data, my_solution, my_col); /*my_solution = 0 vector*/
	dw_printf(IMPORTANCE_DIAG,"Thread %d: my_obj = %3.1f\n", id, my_data->obj);

	/* Multiply D matrix by my_solution and store in vector y. */
	prepare_column(my_solution, y, my_data);

	/* Let master know that I've completed my initial solve. */
	signal_availability(my_data);

	max_iterations = my_data->globals->max_phase1_iterations +
						my_data->globals->max_phase2_iterations;

	/* Now subproblem loops, continually re-solving until signaled by master
	 * to stop or max_iterations is reached.
	 */
	for( j = 1; j < max_iterations; j++ ) {
		/* Throughout this loop, 'col' is used, but in reality, we are only
		 * concerned with col-1 elements.  The 0th element is generally
		 * ignored. */

		pthread_mutex_lock(&sub_data_mutex[id]);

		/* row_duals' * D = (D' * row_duals)' = y  |  q'D = y */
#ifdef USE_INTEL_MKL
		mkl_dcoomv(&yes_trans, &(D->rows), &(my_data->num_cols), &double_one,
				matdescra, my_data->D_vals, my_data->D_row_coords,
				my_data->D_col_coords, &(my_data->D_nnz),
				(my_data->md->row_duals+1), &double_zero, y);
#else
		dw_dcoogemv(my_data->num_cols, D->rows, my_data->D_vals,
				my_data->D_col_coords, my_data->D_row_coords, my_data->D_nnz,
				my_data->md->row_duals+1, y);
#endif

		/* my_data->c ===> d_vector */
#ifdef USE_INTEL_MKL
		if( my_data->phase_one )
			dcopy(&my_col, col_zero_vector, &increment, d_vector, &increment);
		else dcopy(&my_col, my_data->c, &increment, d_vector, &increment);
#else
		if( my_data->phase_one )
			dw_dcopy(my_col, col_zero_vector, d_vector);
		else dw_dcopy(my_col, my_data->c, d_vector);
#endif
		/* Get the elements of y that correspond to my subproblem.
		   Just using my_data->val as a temporary holding place. */
		for( i = 1; i < my_col; i++ ) {
			my_data->double_vector[i] = y[i-1];
		}

		/* -q*D + c = my objective coeffs (stored in d_vector) */
#ifdef USE_INTEL_MKL
		daxpy(&my_col, &double_neg_one, my_data->val, &increment,
				d_vector, &increment);
#else
		dw_daxpy(my_col, double_neg_one, my_data->double_vector, d_vector);
#endif
		for( i = 1; i < my_col; i++ ) {
			glp_set_obj_coef(lp, i, d_vector[i]);
		}

		/* Re-solve based on new objective. */
		dw_printf(IMPORTANCE_DIAG," %d: Going to solve...\n", id);
		ret = glp_simplex(lp, simplex_control_params);
		my_data->obj = glp_get_obj_val(lp);
		if( my_data->globals->enforce_sub_integrality ) {
			glp_intopt(lp, int_parm);
			my_data->obj = glp_mip_obj_val(lp);
		}
		dw_printf(IMPORTANCE_DIAG,
				"Thread %d: Finished solve.  Simplex returns %d. glp_status = %d\n", id, ret, glp_get_status(lp));
		my_data->unbounded = glp_get_status(lp);
		glp_create_index(lp);
		if( my_data->unbounded == GLP_UNBND) {
			//printf("UNBOUNDED: glp_get_unbnd_ray = %d, variable = %s\n",
			//		glp_get_unbnd_ray(lp), glp_get_col_name(lp, glp_get_unbnd_ray(lp)));
			printf("This version of dwsolver does not support unbounded subproblems.\n");
			//my_data->val[glp_get_unbnd_ray(lp)] = 1.0;
			if( glp_get_unbnd_ray(lp) <= glp_get_num_rows(lp)) {
				//my_data->unbounded = GLP_UNBND + 1;
			}
			for( i = 1; i < my_col; i++ ) {
				//printf("********   %s status is %d\n", glp_get_col_name(lp, i), glp_get_col_stat(lp, i));
			}
			for( i = 1; i <= glp_get_num_rows(lp); i++ ) {
				//printf("********   %s status is %d\n", glp_get_row_name(lp, i), glp_get_row_stat(lp, i));
			}
		}

		dw_printf(IMPORTANCE_DIAG,
				"Thread %d: Going to try getting my solution together.\n", id);

		/* Translate solution indices to master indices. */
		dw_printf(IMPORTANCE_DIAG,
				"Thread %d: There are %d cols.\n", id, my_col-1);
		organize_solution(my_data, my_solution, my_col); /* my_solution = 0 */
		dw_printf(IMPORTANCE_DIAG,
				"Thread %d: my_obj = %3.1f\n", id, my_data->obj);

		/* Dx = y  where x is my solution. */
		prepare_column(my_solution, y, my_data);

		/*fflush(stdout); Not a very thread safe function. */
		pthread_mutex_unlock(&sub_data_mutex[id]);
		if( signal_availability(my_data) ) break;
		//sprintf(local_buffer, "sub_%d_%d.cpxlp", id, j);
		//lpx_write_cpxlp(my_data->lp, local_buffer);

	}

	/* Clean up after myself.  my_data will be freed by master. */
	free(y);
	free(d_vector);
	free(simplex_control_params);
	free(my_solution);
	free(local_buffer);
	free(col_zero_vector);
	free(int_parm);
	dw_printf(IMPORTANCE_DIAG,"Thread %d: Exiting.\n", id);
	pthread_exit(NULL);
}

/* This function lets the master know this subproblem is ready for processing on
 * this iteration.  Subproblem puts itself in a queue (service_queue) and then
 * waits ("customers" semaphore).
 */
int signal_availability(subprob_struct* my_data) {
	int rc = 0;
	pthread_mutex_lock(&service_queue_mutex);
	my_data->globals->service_queue[my_data->globals->tail_service_queue] =
		my_data->my_id;
	my_data->globals->tail_service_queue++;
	my_data->globals->tail_service_queue %= my_data->globals->num_clients;
#ifdef USE_NAMED_SEMAPHORES
	sem_post(customers);
#else
	sem_post(&customers);
#endif
	pthread_mutex_unlock(&service_queue_mutex);

	/* Wait for signal here... */
	pthread_mutex_lock(&next_iteration_mutex);
	if (signals->current_iteration == my_data->local_iteration) {
		pthread_cond_wait(&next_iteration_cv, &next_iteration_mutex);
	}
	pthread_mutex_unlock(&next_iteration_mutex);

	/* Should check the my_data->command */
	/* If it's COMMAND_STOP  then the master doesn't need anymore from us? */
	pthread_mutex_lock(&sub_data_mutex[my_data->my_id]);
	if( my_data->command == COMMAND_STOP ) {
		dw_printf(IMPORTANCE_DIAG,
				"%d: Rec'd 'stop' signal.  Bailing right now.  Bye.\n",
				my_data->my_id);
		rc = 1;
	}
	pthread_mutex_unlock(&sub_data_mutex[my_data->my_id]);
	my_data->local_iteration++;
	return rc;
}

/* There used to be a lot more to this function.  Now it simply copies data from
 * one array to another.
 */
void organize_solution(subprob_struct* data, double* sol, int my_col) {
	int i;

	/* Check if master is going to add a convexity variable (lambda). If so,
	 * make sure the convexity constraint is bounded correctly.
	 */
	if( data->r - data->obj > 0.0+TOLERANCE && data->unbounded != GLP_UNBND ) {
		// This would be a good place to set flags or do other things to
		// handle an unbounded ray.

		//dw_printf(IMPORTANCE_DIAG, "Setting row %d bounds to 1.\n",
		//		D->rows + data->my_id + 1);
		//glp_set_row_bnds(master_lp, D->rows + data->my_id + 1, GLP_FX, 1.0, 1.0);
	}
	for( i = 1; i < my_col; i++ ) {
		if( data->globals->enforce_sub_integrality ) {
			sol[i] = glp_mip_col_val(data->lp, i);
		}
		else {
			sol[i] = glp_get_col_prim(data->lp, i);
		}
	}
}

/* This function takes the solution to the subproblem and translates it into a
 * column for the master problem using a pre-computed look-up table.  Also
 * calculates the master problem's objective coeff. for this column.
 */
void prepare_column(double* solution, double* y, subprob_struct* data) {
	int count = 0, i;
	int id = data->my_id;
	pthread_mutex_lock(&glpk_mutex);
	dw_printf(IMPORTANCE_DIAG,"%d: in prepare_column()\n", data->my_id);
	for( i = 1; i <= data->num_cols; i++ ) {
		dw_printf(IMPORTANCE_DIAG," %d:   %s: \t%3.2f\n", id,
				glp_get_col_name(data->lp, i),
				glp_get_col_prim(data->lp, i));
	}

	pthread_mutex_unlock(&glpk_mutex);

	/*   D * solution = y */
#ifdef USE_INTEL_MKL
	mkl_dcoomv(&no_trans, &(D->rows), &(data->num_cols), &double_one, matdescra,
			data->D_vals, data->D_row_coords, data->D_col_coords,&(data->D_nnz),
			solution+1, &double_zero, y);
#else
	dw_dcoogemv(D->rows,  data->num_cols, data->D_vals,
			data->D_row_coords, data->D_col_coords, data->D_nnz, solution+1,
			y);
#endif

	/* Put y into compact form (non-zero coefficients and their indices) */
	for( i = 1; i <= D->rows; i++ ) {
		if( y[i-1] != 0.0 ) {
			count++;
			data->val[count] = y[i-1];
			data->ind[count] = i;
			if( data->globals->verbosity >= OUTPUT_ALL )
				printf("Thread %d: val[%d] = %f, ind[%d] = %d\n",
						data->my_id, count, y[i-1], count, i);
		}
	}
	data->len = count;

	/* Just for debugging purposes, print out some information. */
	if( data->globals->verbosity >= OUTPUT_ALL ) {
		for( i = 0; i < D->rows; i++ )
			printf("Thread %d: y[%d] = %3.1f\n", id, i, y[i]);
		printf("Thread %d: my_data->r = xx, objective = xx\n",
				data->my_id /*data->r, data->obj*/);
		for( i = 1; i <= data->num_cols; i++ ) {
			printf("data->c[%d] = %3.1f, my_solution[%d] = %3.1f (%s)\n", i,
					data->c[i], i, solution[i], glp_get_col_name(data->lp, i));
		}
	}

	/* Subproblem's solution * original cost vector = master objective coeff. */
#ifdef USE_INTEL_MKL
	data->new_obj_coeff = ddot(&(data->num_cols), data->c+1, &increment,
			solution+1, &increment);
#else
	data->new_obj_coeff = dw_ddot(data->num_cols, data->c+1, solution+1);
#endif
	dw_printf(IMPORTANCE_DIAG,
		"Thread %d: my_data->new_obj_coeff = %3.2f\n", id,
		data->new_obj_coeff);
}

/* This is a debugging function that is of dubious use and quality. */
void printD() {
	int i, j, begin, end, count;
	count = 0;
	for( i = 0; i < D->rows; i++ ) {

		begin = D->pointerB[i];
		end   = D->pointerE[i];
		printf("Begin: %d, End %d\n", begin, end);
		for( j = 0; j < (end-begin); j++) {
			printf("%d:%3.1f ", D->columns[count], D->values[count]);
			count++;
		}

		printf("\n");
	}
}
