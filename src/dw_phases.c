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

/*
 * phases.c
 *
 *  Created on: Jun 3, 2009
 *      Author: jrios
 *
 *  The two functions here are nearly identical.  They could/should be a single
 *  function with perhaps an additional parameter passed to it.  The fact that
 *  they work is enough reason right now to NOT change them.
 */
#include <stdlib.h>
#include <stdio.h>
#include <pthread.h>
#include "dw.h"
#include "dw_subprob.h"
#include "dw_support.h"
#include "dw_phases.h"

int phase_1_iteration(subprob_struct* sub_data, faux_globals* fg, int first_run,
		char** obj_names, double* obj_coefs, int* obj_count, master_data* md) {

	int count = 0;
	int index = -1;
	int col;
	int simplex_rc;
	int i;
	static int simplex_iterations = 0;
	static int iteration_count = 0;
	int rc = 0;
	int num_clients = fg->num_clients;
	/* Use these local pointers to grab the y-columns created earlier.  These
	 * will be columns with only one non-zero, indexed by 1.  To be safe,
	 * provide one extra space in memory for this column, thus size of 3. */
	int* ind_local = NULL;
	double* val_local = NULL;
	double* y_accumulators = NULL;
	char* col_name = NULL;
	if( first_run ) {
		ind_local      = (int*)    malloc(sizeof(int)*3);
		val_local      = (double*) malloc(sizeof(double)*3);
		y_accumulators = (double*) calloc(D->rows+1 , sizeof(double));
		col_name       = malloc(sizeof(char)*BUFF_SIZE);
	}

	dw_printf(IMPORTANCE_DIAG, "######## In master_phase_one().\n");

	/* We need a response from each of the clients. 'count' keeps track. */
	while(count < num_clients) {
#ifdef USE_NAMED_SEMAPHORES
		sem_wait(customers);
#else
		sem_wait(&customers);
#endif
		pthread_mutex_lock(&service_queue_mutex);
		index = fg->service_queue[fg->head_service_queue];
		fg->head_service_queue++;
		fg->head_service_queue %= num_clients;
		pthread_mutex_unlock(&service_queue_mutex);
		dw_printf(IMPORTANCE_DIAG, "######## Master servicing thread %d.\n", index);

		count++;

		/* Add column when sub_data[index].r > sub_data[index].obj */
		pthread_mutex_lock(&sub_data_mutex[index]);
		dw_printf(IMPORTANCE_DIAG, "I think r = %3.6f and obj = %3.6f and unbounded = %s\n",
				sub_data[index].r, sub_data[index].obj,
				sub_data[index].unbounded == GLP_UNBND ? "GLP_UNBOUND" : "not GLP_UNBND");

		if( sub_data[index].r - sub_data[index].obj > TOLERANCE || sub_data[index].unbounded == GLP_UNBND) {
			col = glp_add_cols(master_lp, 1);

			if( sub_data[index].unbounded == GLP_UNBND ) {
				printf("Unbounded subproblems are not supported.\n");
				glp_set_col_bnds(master_lp, col, GLP_LO, 0.0, DW_INFINITY);
				glp_set_col_kind(master_lp, col, GLP_CV);
			}
			else {
				glp_set_col_bnds(master_lp, col, GLP_DB, 0.0, 1.0);
				glp_set_col_kind(master_lp, col, GLP_IV);
			}
			obj_names[*obj_count] = malloc(sizeof(char)*BUFF_SIZE);

			/* Lock here is overkill, but it stops warnings from helgrind. */
			pthread_mutex_lock(&next_iteration_mutex);
			if( sub_data[index].unbounded != GLP_UNBND )
				snprintf(obj_names[*obj_count], BUFF_SIZE - 1,
						"lambda_%d_%d", index, signals->current_iteration);
			else
				snprintf(obj_names[*obj_count], BUFF_SIZE - 1,
						"theta_%d_%d", index, signals->current_iteration);
			pthread_mutex_unlock(&next_iteration_mutex);

			glp_set_col_name(master_lp, col, obj_names[*obj_count]);

			dw_printf(IMPORTANCE_DIAG,
					"Master is adding a column called %s.\n",
					obj_names[*obj_count]);

			dw_printf(IMPORTANCE_DIAG,
					"Obj coeff is %f.\n",
					sub_data[index].new_obj_coeff);

			obj_coefs[*obj_count] = sub_data[index].new_obj_coeff;
			*obj_count = *obj_count + 1;

			/* These lines are for the convexity constraint. */
			if( sub_data[index].unbounded != GLP_UNBND ) {
				i = ++(sub_data[index].len); // Add one more row.
				sub_data[index].ind[i] = (D->rows + index + 1); // Index of constraint
				sub_data[index].val[i] = 1.0; // Coefficient of 1.0.
			}
			else {
				i = sub_data[index].len;
				printf("Unbounded subproblem %d is not going to be handled correctly.\n", index);
				fflush(stdout);
			}

			/* Now actually add the new column to constraint matrix. */
			//printf("sub_data[%d].ind[1] = %d\n", index, sub_data[index].ind[1]);
			fflush(stdout);
			if( sub_data[index].ind[1] < 1 ) getc(stdin);
			glp_set_mat_col(master_lp, col, i, sub_data[index].ind, sub_data[index].val);
			glp_set_obj_coef(master_lp, col, 0.0);

			if( first_run )
				for( i = 1; i < sub_data[index].len; i++ ) {
					//printf("Adding %f to y_acc[%d]\n", sub_data[index].val[i], sub_data[index].ind[i]);
					//fflush(stdout);
					y_accumulators[sub_data[index].ind[i]] += sub_data[index].val[i];
					//printf("Adding %f to y_acc[%d]\n", sub_data[index].val[i], sub_data[index].ind[i]);
				}

			fg->x[index][signals->current_iteration] = malloc(sizeof(double)*sub_data[index].num_cols_plus);

			for( i = 1; i < sub_data[index].num_cols_plus; i++ ) {
				fg->x[index][signals->current_iteration][i] =
					sub_data[index].current_solution[i];
			}

			rc++;
		}
		else {
			fg->x[index][signals->current_iteration] = NULL;
		}
		pthread_mutex_unlock(&sub_data_mutex[index]);
		dw_printf(IMPORTANCE_DIAG,
				"######## Master done servicing thread %d.\n", index);

	}
	/* At this point, all subproblems are waiting for my signal to go again. */
	for(count = 0; count < num_clients; count++) {
		/* This doesn't actually make the subproblems go.  Just setting a flag
		 * for when they truly are signaled to continue. */
		sub_data[count].command = COMMAND_GO;
	}

	/* Set the sign of the auxiliary variable as appropriate on the first run
	 * through phase 1.
	 */
	if( first_run ) {
		glp_create_index(master_lp);
		for(i = 1; i < md->rows + 1; i++) {
			if( glp_get_row_lb(master_lp, i) - y_accumulators[i] < 0.0 ) {
				snprintf(col_name, BUFF_SIZE-1, "y_%d", i);
				//printf("Trying to set %s appropriately.\n", col_name);
				col = glp_find_col(master_lp, col_name);
				//printf("Found column to be %d\n", col);
				if( col == 0 ) continue;
				ind_local[1] = i;
				val_local[1] = -1.0;
				//printf("In first_run inner branch.\n");
				//fflush(stdout);
				glp_set_mat_col(master_lp, col, 1, ind_local, val_local);
			}
		}
	}
	/* Solve master again with all of these new columns. */
	simplex_rc = glp_simplex(master_lp, simplex_control_params);
	iteration_count++;
	if( fg->verbosity >= OUTPUT_VERBOSE ) {
		printf("The simplex return code was %d.\n", simplex_rc);
		printf("glp_get_status()    returns %d.\n", glp_get_status(master_lp));
		printf("glp_get_prim_stat() returns %d.\n", glp_get_prim_stat(master_lp));
		printf("glp_get_dual_stat() returns %d.\n", glp_get_dual_stat(master_lp));
	}
	/* Based on the number of simplex iterations and columns added,
	 * deduce if we need to go on, need to stop, or are cycling.
	 */
	if( simplex_iterations < lpx_get_int_parm(master_lp, LPX_K_ITCNT) ) {
		simplex_iterations = lpx_get_int_parm(master_lp, LPX_K_ITCNT);
		if( fg->verbosity >= OUTPUT_VERBOSE)
			printf("There have been %d simplex iterations thus far.\n",
					simplex_iterations);
	}
	/* First run may have no simplex iterations, but other DW iterations should
	 * generate simplex activity.
	 */
	else if(!first_run) {
		/* Maybe good to take some action here? */
		if( fg->verbosity >= OUTPUT_NORMAL ) {
			printf("There weren't any simplex iterations.  I seriously ");
			printf("doubt any further progress.\n");
			//printf("BUT, I'm going to keep going anyway...\n");
		}
	}

	if( fg->verbosity >= OUTPUT_NORMAL ) {
		printf("################################################\n");
		printf("####  Master objective value = %.8e \n", glp_get_obj_val(master_lp));
		printf("################################################\n");
	}

	/* Rebuild the index since we should have new columns with new names. */
	glp_create_index(master_lp);
	if( fg->write_bases ) {
		/* '+1000' to differentiate between phase I and II iterations. */
		write_basis(iteration_count+1000);
	}
	if( glp_get_status(master_lp) == GLP_INFEAS )  {
		/* Need to take some sort of action here? */
		printf("Master is infeasible.\n");
	}
	if( fg->verbosity >= OUTPUT_ALL ) {
		for( i = 1; i <= glp_get_num_rows(master_lp); i++ ) {
			printf("Row %d has dual value %3.1f\n",
					i, glp_get_row_dual(master_lp, i));
		}
	}
	/* Save the dual cost vector for the subprobs' use in next iteration.*/
	for( i = 1; i <= D->rows; i++ ) {
		md->row_duals[i] = glp_get_row_dual(master_lp, i);
	}

	/* This is the 'target' for the subproblem to 'beat' on next iteration.
	 * This target is the dual value of the convexity constraint associated
	 * with the subproblem. */
	for( i = 0; i < num_clients; i ++ ) {
		pthread_mutex_lock(&sub_data_mutex[i]);
		sub_data[i].r = glp_get_row_dual(master_lp, D->rows + 1 + i);
		pthread_mutex_unlock(&sub_data_mutex[i]);
	}

	/* In case a subthread hasn't actually blocked on the next_iteration_cv
	 * signal, we keep track of the current iteration for the subthreads to
	 * check against to see if they are behind and can bypass waiting for the
	 * signal (which they missed). */
	pthread_mutex_lock(&next_iteration_mutex);
	signals->current_iteration++;
	pthread_cond_broadcast(&next_iteration_cv);
	pthread_mutex_unlock(&next_iteration_mutex);

	if( first_run ) {
		free(ind_local);
		free(val_local);
		free(y_accumulators);
		free(col_name);
	}

	dw_printf(IMPORTANCE_DIAG, "######## Leaving phase_1_iteration().\n");
	return rc;
}


/* All of the heavy algorithmic lifting is done here.  A 'master iteration'
 * consists of waiting for each subproblem to offer a column to the master.  As
 * the master receives these columns, it will incorporate the column to the
 * reduced master problem if it is expected to reduce the objective function
 * obtained from the previous iteration.  If the column is included, the actual
 * values of the variables associated with the subproblem are stored so that
 * they can be easily retrieved later when forming the final solution.  After
 * all subproblems have reported, the reduced master is solved to optimality.
 * This solving process will be a number of simplex iterations equal to the
 * number of columns added to the problem.  The master may or may not be coded
 * to enforce integrality of solutions at this point.  The master also signals
 * all subproblems to go solve again after setting the row dual values (these
 * are necessary for the subproblems to set their new objective functions.
 */
int phase_2_iteration(subprob_struct* sub_data, faux_globals* fg, master_data* md) {
	int count = 0;
	int index = -1;
	int col;
	int i;
	int rc = 0;
	static int iteration_count = 0;
	static int simplex_iterations = 0;
	int j;
	int simplex_rc;
	int num_clients = fg->num_clients;
	char* buffer = (char*) malloc(sizeof(char)*BUFF_SIZE);
	//glp_iocp* int_parm = malloc(sizeof(glp_iocp));
	//glp_init_iocp(int_parm);

	dw_printf(IMPORTANCE_DIAG, "######## In master_iteration().\n");
	/* We need a response from each of the clients. 'count' keeps track. */
	while(count < num_clients) {
#ifdef USE_NAMED_SEMAPHORES
		sem_wait(customers);
#else
		sem_wait(&customers);
#endif
		pthread_mutex_lock(&service_queue_mutex);
		index = fg->service_queue[fg->head_service_queue];
		fg->head_service_queue++;
		fg->head_service_queue %= num_clients;
		pthread_mutex_unlock(&service_queue_mutex);
		dw_printf(IMPORTANCE_DIAG, "######## Master servicing thread %d.\n", index);
		count++;

		/* Add column when sub_data[index].r > sub_data[index].obj */
		pthread_mutex_lock(&sub_data_mutex[index]);
		dw_printf(IMPORTANCE_DIAG, "I think r = %3.6f and obj = %3.6f and unbounded = %s\n",
				sub_data[index].r, sub_data[index].obj,
				sub_data[index].unbounded == GLP_UNBND ? "GLP_UNBOUND" : "not GLP_UNBND");

		if( sub_data[index].r - sub_data[index].obj > 0.0 + TOLERANCE || sub_data[index].unbounded == GLP_UNBND) {
			col = glp_add_cols(master_lp, 1);

			if( sub_data[index].unbounded == GLP_UNBND ) {
				printf("Unbounded subproblems are not supported.\n");
				glp_set_col_bnds(master_lp, col, GLP_LO, 0.0, DW_INFINITY);
				glp_set_col_kind(master_lp, col, GLP_CV);
			}
			else {
				glp_set_col_bnds(master_lp, col, GLP_DB, 0.0, 1.0);
				glp_set_col_kind(master_lp, col, GLP_IV);
			}
			pthread_mutex_lock(&next_iteration_mutex);
			if( sub_data[index].unbounded != GLP_UNBND )
				snprintf(buffer, BUFF_SIZE - 1,
						"lambda_%d_%d", index, signals->current_iteration);
			else
				snprintf(buffer, BUFF_SIZE - 1,
						"theta_%d_%d", index, signals->current_iteration);
			pthread_mutex_unlock(&next_iteration_mutex);
			dw_printf(IMPORTANCE_DIAG, "Master is adding a column called %s.\n",
					buffer);

			if( sub_data[index].unbounded != GLP_UNBND ) {
				/* These next 3 lines are for the convexity constraint. */
				i = ++(sub_data[index].len); // Add one more row.
				sub_data[index].ind[i] = (D->rows + index + 1); // Index of constraint
				sub_data[index].val[i] = 1.0; // Coeff of 1.0.
			}
			else
				i = sub_data[index].len;

			/* Now actually add the new column to constraint matrix. */
			glp_set_mat_col(master_lp, col, i, sub_data[index].ind, sub_data[index].val);
			glp_set_obj_coef(master_lp, col, sub_data[index].new_obj_coeff);
			glp_set_col_name(master_lp, col, buffer);

			if( glp_get_col_prim(master_lp, col) < (1.0-TOLERANCE) && glp_get_col_prim(master_lp, col) > TOLERANCE) {
				printf("!!!!!  %s was just added and its value is %3.5f\n", buffer, glp_get_col_prim(master_lp, col) );
				//Get column.  Step through the indices and check aux var value associated with row.
				int* ind = malloc(sizeof(int)*glp_get_num_rows(master_lp));

				double* val = NULL;
				int len = glp_get_mat_col(master_lp, col, ind, val);
				for( j = 1; j <= len; j++) {
					if( glp_get_row_prim(master_lp, ind[j]) <= TOLERANCE )  {
						int* ind2 = malloc(sizeof(int)*glp_get_num_cols(master_lp));
						double* val2 = malloc(sizeof(double)*glp_get_num_cols(master_lp));
						printf("      %s is at its bound (%2.1f).\n", glp_get_row_name(master_lp, ind[j]), glp_get_row_ub(master_lp, ind[j]));
						int len2 = glp_get_mat_row(master_lp, ind[j], ind2, val2);

						int k;
						for( k = 1; k <= len2; k++ ) {
							if( glp_get_col_prim(master_lp, ind2[k]) > TOLERANCE ) {
								printf("          %s = (%2.1f)*%f\n", glp_get_col_name(master_lp, ind2[k]), val2[k], glp_get_col_prim(master_lp, ind2[k]));
							}
						}
						free(ind2); free(val2);
					}
				}
				free(ind);
			}
			rc++;

			fg->x[index][signals->current_iteration] = malloc(sizeof(double)*sub_data[index].num_cols_plus);

			for( i = 1; i < sub_data[index].num_cols_plus; i++ ) {
				fg->x[index][signals->current_iteration][i] =
					sub_data[index].current_solution[i];
			}
		}
		else {
			fg->x[index][signals->current_iteration] = NULL;
		}
		pthread_mutex_unlock(&sub_data_mutex[index]);

#ifdef VERBOSE
		printf("######## Master done servicing thread %d.\n", index);
#endif
	}
	/* At this point, all subproblems are waiting for my signal to go again. */
	for(count = 0; count < num_clients; count++) {
		sub_data[count].command = COMMAND_GO; // This does nothing right now?
	}

	/* Solve master again with all of these new columns. */
	simplex_rc = glp_simplex(master_lp, simplex_control_params);
	iteration_count++;
	if( fg->verbosity >= OUTPUT_VERBOSE ) {
		printf("The simplex return code was %d.\n", simplex_rc);
		printf("glp_get_status()    returns %d.\n", glp_get_status(master_lp));
		printf("glp_get_prim_stat() returns %d.\n", glp_get_prim_stat(master_lp));
		printf("glp_get_dual_stat() returns %d.\n", glp_get_dual_stat(master_lp));
	}

	/* Based on the number of simplex iterations and columns added,
	 * deduce if we need to go on, need to stop, or are cycling.
	 */
	if( simplex_iterations < lpx_get_int_parm(master_lp, LPX_K_ITCNT) ) {
		simplex_iterations = lpx_get_int_parm(master_lp, LPX_K_ITCNT);
		if( fg->verbosity >= OUTPUT_VERBOSE )
			printf("There have been %d simplex iterations thus far.\n", simplex_iterations);
		if( rc <= 0 ) {
			if( fg->verbosity >= OUTPUT_NORMAL )
				printf("No columns added, but simplex made progress.\n");
			rc = -1; /* if rc == 0, then outer loop will bail. */
			//getchar();
		}
	}
	else if ( rc > 0 ){ /* Added columns, but no simplex iterations. */
		rc = iteration_count == 1 ? 1 : 0; /* Force outer loop to bail. */
		if( rc == 0 && fg->verbosity >= OUTPUT_NORMAL) {
			printf("There was no simplex progress made despite columns being ");
			printf("added. Finishing now.\n");
		}
		//rc = 1;
	}
	else { /* Converged. */
		if( fg->verbosity >= OUTPUT_NORMAL )
			printf("I think we've converged on an optimum.\n");
	}

	if( check_col_integrality() ) {

		dw_printf(IMPORTANCE_DIAG, "There are basic variables with non-integer values.\n");
		dw_printf(IMPORTANCE_DIAG, "But I'm not doing anything about it.\n");
//		printf("I want to see if there is an equivalent integer solution.  Solving...\n");
//		simplex_rc = glp_intopt(master_lp, parm);
//		printf("Integer optimization returned %d.\n", simplex_rc);
//		printf("Integer optimal value is %3.2f.\n", glp_mip_obj_val(master_lp));
	}
	else dw_printf(IMPORTANCE_DIAG, "The basic variables are all integer.\n");

	//if( iteration_count % 5 == 0 ) purge_nonbasics();

	//purge_nonbasics();

	//printf("Checking aux vars for binding constraints...\n");
	//check_aux_vars();

	/* Rebuild the index since we should have new columns with new names. */
	glp_create_index(master_lp);
	if( fg->write_bases ) write_basis(iteration_count);
	if( glp_get_status(master_lp) == GLP_INFEAS ) printf("Master is infeasible.\n");
	for( i = 1; i <= glp_get_num_rows(master_lp); i++ ) {
		//printf("Row %d has dual value %3.1f\n", i, glp_get_row_dual(master_lp, i));
	}

	for( i = 1; i <= glp_get_num_cols(master_lp); i++ ) {
		if( glp_get_col_stat(master_lp, i) == GLP_BS ) {
			//printf("Col %d (%s) is basic with optimal value %3.2f and obj. coeff of %3.2f\n",
			//		i, glp_get_col_name(master_lp, i), glp_get_col_prim(master_lp, i), glp_get_obj_coef(master_lp, i));
		}
	}


	if( fg->verbosity >= OUTPUT_NORMAL ) {
		printf("################################################\n");
		printf("####  Master objective value = %e \n", glp_get_obj_val(master_lp));
		printf("################################################\n");
	}
	/* Save the dual cost vector for the subprobs' use in next iteration.*/
	pthread_mutex_lock(&master_mutex);
	for( i = 1; i <= D->rows; i++ ) {
		md->row_duals[i] = glp_get_row_dual(master_lp, i);
	}
	pthread_mutex_unlock(&master_mutex);

	/* This is the target for the subproblem to beat on next iteration.
	 * This target is the dual value of the convexity constraint associated
	 * with the subproblem. */
	for( i = 0; i < num_clients; i ++ ) {
		pthread_mutex_lock(&sub_data_mutex[i]);
		sub_data[i].r = glp_get_row_dual(master_lp, D->rows + 1 + i);
		pthread_mutex_unlock(&sub_data_mutex[i]);
	}

	/* In case a subthread hasn't actually blocked on the next_iteration_cv
	 * signal, we keep track of the current iteration for the subthreads to
	 * check against to see if they are behind and can bypass waiting for the
	 * signal (which they missed). */
	pthread_mutex_lock(&next_iteration_mutex);
	signals->current_iteration++;
	pthread_cond_broadcast(&next_iteration_cv);
	pthread_mutex_unlock(&next_iteration_mutex);

	free(buffer);

	dw_printf(IMPORTANCE_DIAG, "######## Leaving phase_2_iteration().\n");
	return rc;
}
