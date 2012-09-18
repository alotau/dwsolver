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
 * dw_rounding.c
 *
 *  Created on: Dec 15, 2008
 *      Author: jrios
 */


#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <pthread.h>
#include <math.h>

#ifdef USE_INTEL_MKL
#include <mkl_spblas.h>  /* Sparse matrix stuff. */
#include <mkl.h>         /* Other Math Kernel Library stuff. */
#endif

#include "dw_blas.h"
#include "dw.h"
#include "dw_subprob.h"
#include "dw_support.h"
#include "dw_rounding.h"
#include "dw_phases.h"

static int post_round_zeros   = 0;
static int post_round_ones    = 0;
static int pre_round_zeros    = 0;
static int pre_round_ones     = 0;
static int pre_round_nonzeros = 0;

static int round_up = 1;
static double round_param = 0.5;

int process_solution(subprob_struct* sub_data, char* out_file_name, int mode) {
	int i;
	int j;
	int rc;
	int subprob;
	int iteration;
	int num_clients = sub_data->globals->num_clients;
	int* num_infos = malloc(sizeof(int)*num_clients);
	const char* col_name;
	char* local_col_name = malloc(sizeof(char)*BUFF_SIZE);
	double var_value;
	void *status;
	FILE* out_file;

	/* For convenience. */
	faux_globals* globals = sub_data->globals;

	sol_info* temp_info;
	sol_info** solution_organization = malloc(sizeof(sol_info*)*num_clients);
	int_thread_data* t_data = malloc(sizeof(int_thread_data)*num_clients);
	pthread_t* threads = (pthread_t*) malloc(sizeof(pthread_t)*num_clients);

	globals->final_x   =
		calloc(glp_get_num_cols(original_master_lp)+1, sizeof(double));
	globals->relaxed_x =
		calloc(glp_get_num_cols(original_master_lp)+1, sizeof(double));

	dw_printf(IMPORTANCE_DIAG, "Entered round_solution()...\n");

	/* Initialize some data structures. */
	for( i = 0; i < num_clients; i++ ) {
		num_infos[i] = 0;
		solution_organization[i] = NULL;
		sub_data[i].globals->final_x = sub_data->globals->final_x;
		sub_data[i].globals->relaxed_x = sub_data->globals->relaxed_x;
	}

	/* Examine each convexity constraint. Locate the subproblem and iteration
	 * associated with each non-zero valued variable in the convexity row.
	 */
	dw_printf(IMPORTANCE_DIAG, "Going into row for loop...\n");

	for( i = D->rows; i <= glp_get_num_cols(master_lp); i++ ) {
		var_value = glp_get_col_prim(master_lp, i);
		dw_printf(IMPORTANCE_DIAG, "Col %s has primal value %f\n", glp_get_col_name(master_lp, i), var_value);
		if( var_value != 0.0 ) {

			/* It will help to recall that the convexity constraints are named
			 * for the subproblem and iteration: lambda_subprob_iteration.
			 * It is because of this naming convention that we parse as follows:
			 */
			col_name = glp_get_col_name(master_lp, i);
			if( strncmp("lambda", col_name, strlen("lambda")) != 0 ) continue;
			//dw_printf(IMPORTANCE_DIAG, "Col = %d (%s).\n", i, col_name);
			strcpy(local_col_name, col_name);
			/* Tease out the subprob number and the iteration. */
			strtok(local_col_name, "_");
			/* Would probably make sense to do safety checks in here. */
			subprob   = atoi(strtok(NULL, "_"));
			iteration = atoi(strtok(NULL, "_"));
			/* Have we encountered a lambda from this subproblem before? */
			if( num_infos[subprob] == 0 ) {
				/* Begin linked list of non-zero vars for this subproblem. */
				solution_organization[subprob] = malloc(sizeof(sol_info));
				temp_info = solution_organization[subprob];
			}
			else {
				/* Add to linked list of non-zero vars for this subproblem. */
				temp_info = solution_organization[subprob];
				for( j = 0; j < num_infos[subprob]-1; j++ ) {
					temp_info = temp_info->next;
				}
				temp_info->next = malloc(sizeof(sol_info));
				temp_info = temp_info->next;
			}
			/* Now populate the sol_info data structure. */
			temp_info->iteration = iteration;
			temp_info->subprob = subprob;
			temp_info->value = var_value;
			temp_info->next = NULL;

			num_infos[subprob]++;
		}
	}

	/* There is probably a better way to handle this... */
	if( (out_file = fopen(out_file_name, "w")) == NULL ) {
		dw_printf(IMPORTANCE_VITAL,"Problem opening file: %s\n", out_file_name);
		dw_printf(IMPORTANCE_VITAL,"Aggressively bailing...\n");
		return -1;
	}

	/* Now launch threads to handle the rounding calculations.  One thread per
	 * subproblem. */
	dw_printf(IMPORTANCE_DIAG, "Going into client for loop...\n");
	for( i = 0; i < num_clients; i++ ) {
		//sub_data->globals->x[i][signals->current_iteration] =
		//	malloc(sizeof(double)*sub_data[i].num_cols_plus);
		if( num_infos[i] == 0 ) solution_organization[i] = NULL;
		t_data[i].si = solution_organization[i];
		t_data[i].num_si = num_infos[i];
		t_data[i].id = i;
		t_data[i].sub_data = sub_data[i];
		t_data[i].zero_file = out_file;
		/* Actually create the thread now. */
		if( mode == DW_MODE_ROUND_SOL || mode == DW_MODE_INT_SOL ) {
			rc = pthread_create(&threads[i],
					&attr,
					rounding_thread,
					(void *)&t_data[i]);
		}
		else if( mode == DW_MODE_PRINT_RELAXED ) {
			rc = pthread_create(&threads[i],
					&attr,
					solution_printing_thread,
					(void *)&t_data[i]);
		}
		else {
			dw_printf(IMPORTANCE_VITAL, "UNRECOGNIZED MODE!\n");
			exit(-1);
		}
		if (rc) {
			printf("ERROR; return code from pthread_create() is %d\n", rc);
			exit(-1);
		}
	}

	/* This is where we wait for subthreads to join. */
	dw_printf(IMPORTANCE_DIAG, "%-50s\n", "### Going to wait for solution processing subthreads...");
	for(i=0; i<num_clients; i++) {
		rc = pthread_join(threads[i], &status);
		if (rc) {
			printf("ERROR; return code from pthread_join() is %d\n", rc);
			printf("Going to continue waiting for all joins, but behavior now unpredictable.\n");
			//exit(-1);
		}
		if( sub_data->globals->verbosity >= OUTPUT_ALL ) {
			printf("### Completed join with thread %d status = %ld\n",i,
					(long)status);
		}
	}
	dw_printf(IMPORTANCE_DIAG, "DONE!\n");

	/* For sanity, check the number of zero-valued, one-valued, non-integral
	 * variables before and after the rounding procedure.  Note that the
	 * relaxed solution found by the DW algorithm is not modified by the
	 * rounding procedure.  "final_x" holds the result of the rounding.
	 */
	if( mode == DW_MODE_ROUND_SOL || mode == DW_MODE_INT_SOL ) {
		if( sub_data->globals->verbosity >= OUTPUT_NORMAL ) {
			printf("# of pre-rounding zeros was %d\n", pre_round_zeros);
			printf("# of pre-rounding ones was  %d\n", pre_round_ones);
			printf("# of pre-rounding non-ints was %d\n", pre_round_nonzeros);
			//	}
			check_broken_constraints(sub_data[0].globals->relaxed_x);

			//	if( sub_data->globals->verbosity >= OUTPUT_NORMAL ) {
			printf("# of post-rounding zeros is %d\n",
					post_round_zeros+pre_round_zeros);
			printf("# of post-rounding ones is %d\n",
					post_round_ones+pre_round_ones);
			check_broken_constraints(sub_data[0].globals->final_x);
		}
	}

	fclose(out_file);
	//free(local_col_name);



	for( i = 0; i < num_clients; i++ ) {
		sol_info* soln_info_ptr_a = solution_organization[i];
		sol_info* soln_info_ptr_b;
		while( soln_info_ptr_a != NULL ) {
			soln_info_ptr_b = soln_info_ptr_a->next;
			free(soln_info_ptr_a);
			soln_info_ptr_a = soln_info_ptr_b;
		}
//		for( j = 0; j < num_infos[i]; j++ ) {
////			sol_info* soln_info_ptr = solution_organization[i];
////			while( soln_info_ptr != NULL && soln_info_ptr->next != NULL )
////				soln_info_ptr = soln_info_ptr->next;
////			if( soln_info_ptr != NULL ) free(soln_info_ptr);
//		}
//
//		free(sub_data->globals->x[i][signals->current_iteration]);
//		free(solution_organization[i]);
	}
	free(threads);
	free(t_data);
	free(num_infos);
	free(local_col_name);
	free(sub_data[0].globals->final_x);
	free(sub_data[0].globals->relaxed_x);
	free(solution_organization);
	return 0;
}

void check_broken_constraints(double* xx) {
	//if( globals->dw_verbosity >= OUTPUT_NORMAL )
		//printf("About to check for broken constraints...\n");
	int i;
	int j;
	int num_over_capacity = 0;
	int*    ind = malloc(sizeof(int)*D->cols+1);
	double* val = malloc(sizeof(double)*D->cols+1);
	val[0] = 0.0;
	ind[0] = 0;
	int nz;
	double product = 0.0;
	double coeff = 0.0;
	for( i = 1; i <= D->rows; i++ ) {
		/* Get the constraint row. */
		nz = glp_get_mat_row(original_master_lp, i, ind, val);

		/* Multiply the row by the values of the variables. */
#ifdef USE_INTEL_MKL
		product = ddoti(&nz, val+1, ind+1, xx+1);
#else
		product = dw_ddoti(nz, val+1, ind+1, xx);
#endif
		if( (product-TOLERANCE) > glp_get_row_ub(original_master_lp, i)) {
			printf(" product for row %d (%s) is %3.2f, but bound is %3.2f\n",
					i, glp_get_row_name(original_master_lp, i), product,
					glp_get_row_ub(original_master_lp, i));
			num_over_capacity++;
			for( j = nz; j > 0; j-- ) {
				//printf("    coeff on %s = %3.1f\n", glp_get_col_name(original_master_lp, ind[j]), val[j]);
				//printf("    that variable has value %3.1f\n", xx[ind[j]]);
			}

		}

	}
	//if( globals->dw_verbosity >= OUTPUT_NORMAL )
	printf("There were %d constraint(s) that were violated.\n",
			num_over_capacity);
	nz = 0;
	for( i = 1; i <= glp_get_num_cols(original_master_lp); i++ ) {
		coeff = glp_get_obj_coef(original_master_lp, i) ;
		if( coeff != 0.0 ) {
			nz++;
			ind[nz] = i;
			val[nz] = coeff;
		}
	}
#ifdef USE_INTEL_MKL
		product = ddoti(&nz, val+1, ind+1, xx+1);
#else
		product = dw_ddoti(nz, val+1, ind+1, xx);
#endif
	//if( globals->dw_verbosity >= OUTPUT_NORMAL )
	printf("The objective is now valued at %3.6f\n",
			product+glp_get_obj_coef(master_lp, 0));
	free(ind);
	free(val);
}

void print_zeros_simple(int_thread_data* my_data) {
	int j;
	int subprob = my_data->id;
	int curr_iter = signals->current_iteration;
	//int time;
	//const char* var_name;
	char* local_col_name = malloc(sizeof(char)*BUFF_SIZE);
	char* curr_flight = malloc(sizeof(char)*BUFF_SIZE);
	char* curr_sector = malloc(sizeof(char)*BUFF_SIZE);


	//printf("in print_zeros_simple...\n");
	strcpy(curr_flight, "COOTIE1");
	strcpy(curr_sector, "COOTIE2");
	/* Print zeros" */
	//printf("%d: Printing zeros...\n", my_data->id);
	for( j = 1; j < my_data->sub_data.num_cols_plus; j++ ) {
		parse_zero_var(my_data->sub_data.globals->x[subprob][curr_iter][j],
				j, my_data->sub_data.lp, my_data->zero_file);
//		if( my_data->sub_data.globals->x[subprob][curr_iter][j] < TOLERANCE) {
//			var_name = glp_get_col_name(my_data->sub_data.lp, j);
//			strcpy(local_col_name, var_name);
//			if( strtok(local_col_name, "(") == NULL ) printf("NULL 1");
//			strcpy( curr_flight, strtok(NULL, ",") );
//			strcpy( curr_sector, strtok(NULL, ",") );
//			time = atoi(strtok(NULL, ")"));
//			fprintf(zero_file, "%s %s\n", curr_flight, curr_sector);
//		}

	}


	//printf("leaving print_zeros_simple...\n");
	free(local_col_name);
	free(curr_flight);
	free(curr_sector);

}

void print_zeros(int_thread_data* my_data) {
	int j;
	int subprob = my_data->id;
	int curr_iter = signals->current_iteration;

	int ground_delay = 0;
	int time;
	int prev_delay  = -1;
	int curr_delay  = -2;
	int new_delay   = -3;

	int num_zeros = 0;

	const char* var_name;
	char* local_col_name = malloc(sizeof(char)*BUFF_SIZE);
	char* curr_flight = malloc(sizeof(char)*BUFF_SIZE);
	char* temp_flight = malloc(sizeof(char)*BUFF_SIZE);
	char* prev_flight = malloc(sizeof(char)*BUFF_SIZE);
	char* curr_sector = malloc(sizeof(char)*BUFF_SIZE);
	char* prev_sector = malloc(sizeof(char)*BUFF_SIZE);
	char* old_sector  = malloc(sizeof(char)*BUFF_SIZE);

	//if( globals->dw_verbosity >= OUTPUT_NORMAL )
		//printf("in print_zeros...\n");
	strcpy(curr_flight, "COOTIE1");
	strcpy(curr_sector, "COOTIE2");
	strcpy(prev_flight, "COOTIE3");
	strcpy(prev_sector, "COOTIE4");
	strcpy(old_sector,  "COOTIE5");
	/* Print zeros" */
	//printf("%d: Printing zeros...\n", my_data->id);
	for( j = 1; j < my_data->sub_data.num_cols_plus; j++ ) {
		if( my_data->sub_data.globals->x[subprob][curr_iter][j] < TOLERANCE) {
			num_zeros++;
			var_name = glp_get_col_name(my_data->sub_data.lp, j);
			if( var_name == NULL ) printf("Returned null (subprob %d, col %d), continuing...\n", subprob, j);
			strcpy(local_col_name, var_name);
			if( strtok(local_col_name, "(") == NULL ) printf("NULL 1");
			strcpy( curr_flight, strtok(NULL, ",") );
			strcpy( curr_sector, strtok(NULL, ",") );
			time = atoi(strtok(NULL, ")"));
			//printf("%s ( %s, %s ) = %3.2f\n", var_name, curr_flight, curr_sector, x[subprob][curr_iter][j]);


			/* NEW FLIGHT */
			if( strcmp(curr_flight, prev_flight) != 0 ) {
				//printf("new flight.\n");
				strcpy(prev_flight, curr_flight);
				ground_delay = 0;
				curr_delay = 1;
				/* Set prev_delay appropriately by looking at previous variable. */
				if( j != 1 ) { /* Not first variable in subproblem? Check previous variable. */
					var_name = glp_get_col_name(my_data->sub_data.lp, j-1);
					strcpy(local_col_name, var_name);
					if( strtok(local_col_name, "(") == NULL ) printf("NULL 2");
					strcpy( temp_flight, strtok(NULL, ",") );
					strcpy( prev_sector, strtok(NULL, ",") );
				}
				else { /* First variable in subproblem.  No previous variable. */
					//printf("here.\n");
					strcpy( temp_flight, "DEFINITELY NOT A MATCH");
					ground_delay = 1;
				}

				if( strcmp(curr_flight, temp_flight) != 0 ) { /* Previous variable refers to different flight. */
					//prev_delay = DW_INFINITY;
					prev_delay = 0;
					ground_delay = 1;
					strcpy(old_sector, "COOTIECOOCOO");
					strcpy(prev_sector, curr_sector);
					//printf("Setting prev_delay to DW_INFINITY.\n");
				}
				else { /* Previous variable referred to same flight (but must have had value of 1.0). */
					prev_delay = 0;
					strcpy(old_sector, prev_sector);
				}
			}


			/* SAME FLIGHT AS BEFORE. */
			else { /* Same flight. */
				//printf("same flight.\n");
				if( strcmp(curr_sector, prev_sector) != 0 ) { /* New sector. */
					if( curr_delay > prev_delay ) { /* New delay. */
						new_delay  = curr_delay - prev_delay;
						prev_delay = curr_delay;
						//printf("NEW DELAY: Flight %s, Sector %s, %d minutes.\n", curr_flight, old_sector, new_delay);
						//printf("%s %s %d\n", curr_flight, old_sector, new_delay);
						//printf("%s %s %d\n", curr_flight, old_sector, new_delay);
						fprintf(my_data->zero_file,"%s %s %d\n", curr_flight, old_sector, new_delay);
					}
					else { /* No new delay. */
						prev_delay = curr_delay;
					}

					if( ground_delay ) {
						ground_delay = 0;
						//printf("NEW_DELAY: Flight %s, Sector %s, %d minutes.\n", curr_flight, prev_sector, curr_delay);
						//printf("%s %s %d\n", curr_flight, prev_sector, curr_delay);
						fprintf(my_data->zero_file,"%s %s %d\n", curr_flight, prev_sector, curr_delay);
					}

					strcpy(old_sector,  prev_sector);
					strcpy(prev_sector, curr_sector);
					curr_delay = 1;
				}
				else { /* Same sector as before. */
					//printf("Incrementing curr_delay.\n");
					curr_delay++;
				}
			}
		}
	}
	//if( globals->dw_verbosity >= OUTPUT_NORMAL )
		//printf("leaving print_zeros...\n");
	fflush(stdout);
	if( num_zeros > 0 )
		printf("num_zeros %d\n", num_zeros);
	free(local_col_name);
	free(curr_flight);
	free(temp_flight);
	free(prev_flight);
	free(curr_sector);
	free(prev_sector);
	free(old_sector);
}

void* solution_printing_thread(void* arg) {
	int i, j, subprob, curr_iter, iteration;
	double val;
	int_thread_data* my_data = (int_thread_data*) arg;
	sol_info* temp_info;

	temp_info = my_data->si;
	subprob = my_data->id;
	curr_iter = signals->current_iteration;

	/* Now use my_data to fill in x[i][current_iteration] with final var values. */
	my_data->sub_data.globals->x[subprob][curr_iter] =
		calloc((my_data->sub_data).num_cols_plus, sizeof(double));

	for(i = 0; i < my_data->num_si; i++ ) {

		iteration = temp_info->iteration;
		val = temp_info->value;

		for( j = 1; j < (my_data->sub_data).num_cols_plus; j++ ) {
			my_data->sub_data.globals->x[subprob][curr_iter][j] +=
				my_data->sub_data.globals->x[subprob][iteration][j] * val;
		}
		temp_info = temp_info->next;
	}
	pthread_mutex_lock(&glpk_mutex);
	for( j = 1; j < (my_data->sub_data).num_cols_plus; j++ ) {
		fprintf(my_data->zero_file, "%s\t%f\n",
				glp_get_col_name(my_data->sub_data.lp, j),
				my_data->sub_data.globals->x[subprob][curr_iter][j]);
	}
	pthread_mutex_unlock(&glpk_mutex);

	free(my_data->sub_data.globals->x[subprob][curr_iter]);
}

void* rounding_thread(void* arg) {

	int i, j;
	int iteration;
	int curr_iter;
	int subprob;

	double val;
	const char* local_col_name = malloc(sizeof(char)*BUFF_SIZE);
	char* curr_flight = malloc(sizeof(char)*BUFF_SIZE);
	char* prev_flight = malloc(sizeof(char)*BUFF_SIZE);
	char* curr_sector = malloc(sizeof(char)*BUFF_SIZE);
	char* prev_sector = malloc(sizeof(char)*BUFF_SIZE);

	strcpy(curr_flight, "COOTIE");
	strcpy(curr_sector, "COOTIE");
	strcpy(prev_flight, "COOTIE");
	strcpy(prev_sector, "COOTIE");

	sol_info* temp_info;
	int_thread_data* my_data = (int_thread_data*) arg;
	pthread_mutex_lock(&master_lp_ready_mutex);
	int my_round = round_up;
	round_up = 1 - round_up;
	//printf("Subproblem %d claiming round_up with value %d\n", my_data->id, my_round);
	pthread_mutex_unlock(&master_lp_ready_mutex);

	double temp;

	subprob = my_data->id;
	curr_iter = signals->current_iteration;

	//printf("Hello from rounding_thread %d\n", my_data->id);
	//printf("There are %d rows in my subproblem and %d variables.\n", glp_get_num_rows(my_data->sub_data.lp), glp_get_num_cols(my_data->sub_data.lp));

	/* Now use my_data to fill in x[i][current_iteration] with final var values. */
	my_data->sub_data.globals->x[subprob][curr_iter] =
		calloc((my_data->sub_data).num_cols_plus, sizeof(double));

	/* Get weighted sum of the appropriate solutions. Store in x[curr_iter]. */
	//printf("%d entering for loop.  Will loop %d times.\n", my_data->id, my_data->num_si);
	temp_info = my_data->si;
	for(i = 0; i < my_data->num_si; i++ ) {

		iteration = temp_info->iteration;
		val = temp_info->value;

		for( j = 1; j < (my_data->sub_data).num_cols_plus; j++ ) {
			my_data->sub_data.globals->x[subprob][curr_iter][j] += my_data->sub_data.globals->x[subprob][iteration][j] * val;
		}
		temp_info = temp_info->next;
	}

	/* Now somehow round them off? */
	for( j = 1; j < (my_data->sub_data).num_cols_plus; j++ ) {
		//if( x[subprob][curr_iter][j] > (0.0+TOLERANCE) ) {
		if( my_data->sub_data.globals->x[subprob][curr_iter][j] > (0.0+TOLERANCE) &&
				my_data->sub_data.globals->x[subprob][curr_iter][j] < (1.0-TOLERANCE) ) {
			/* This var isn't one or zero and, thus, needs integerization. */

			local_col_name = glp_get_col_name(my_data->sub_data.lp, j);
			i = glp_find_col(original_master_lp, local_col_name );

			pthread_mutex_lock(&glpk_mutex);
			if( /*x[subprob][curr_iter][j] == 0.0*/ 0 ) pre_round_zeros++;
			else if(my_data->sub_data.globals->x[subprob][curr_iter][j] == 1.0 )
				pre_round_ones++;
			else if(my_data->sub_data.globals->x[subprob][curr_iter][j] <= 1.0 )
				pre_round_nonzeros++;
			else printf("UH OH! ODD VALUED VARIABLE?  x[%d][%d][%d] = %e\n",
					subprob, curr_iter, j ,
					my_data->sub_data.globals->final_x[i]);
			pthread_mutex_unlock(&glpk_mutex);

			my_data->sub_data.globals->relaxed_x[i] = my_data->sub_data.globals->x[subprob][curr_iter][j];

			temp = my_data->sub_data.globals->x[subprob][curr_iter][j];
			if( my_data->sub_data.globals->x[subprob][curr_iter][j] < round_param )
				my_data->sub_data.globals->x[subprob][curr_iter][j] = 0.0;
			else if(my_data->sub_data.globals->x[subprob][curr_iter][j] >= (1.0 - round_param))
				my_data->sub_data.globals->x[subprob][curr_iter][j] = 1.0;
			else my_data->sub_data.globals->x[subprob][curr_iter][j] = my_round ? 1.0 : 0.0;
			//x[subprob][curr_iter][j] = floor(x[subprob][curr_iter][j] + 0.5);
			//if( temp != x[subprob][curr_iter][j] ) printf("x[%d][%d][%d] was %e and is now %e (my_round = %d)\n", subprob, curr_iter, j ,temp, x[subprob][curr_iter][j], my_round);



			//printf("%s: col %d in sub is %d in master.\n", local_col_name, j, i);
			my_data->sub_data.globals->final_x[i] = my_data->sub_data.globals->x[subprob][curr_iter][j];
			/* Signal all subproblems that master is set up and ready. */
			pthread_mutex_lock(&master_lp_ready_mutex);
			if( my_data->sub_data.globals->final_x[i] == 0.0) post_round_zeros++;
			else if( my_data->sub_data.globals->final_x[i] == 1.0 ) post_round_ones++;
			else printf("UH OH! ROUNDING DIDN'T WORK?  final_x[%d] = %f\n", i, my_data->sub_data.globals->final_x[i]);
			pthread_cond_broadcast(&master_lp_ready_cv);
			pthread_mutex_unlock(&master_lp_ready_mutex);


		//	printf("%e\n", x[subprob][curr_iter][j]);
		}
		else if( my_data->sub_data.globals->x[subprob][curr_iter][j] <= (0.0+TOLERANCE) ) {
			pthread_mutex_lock(&glpk_mutex);
			pre_round_zeros++;
			pthread_mutex_unlock(&glpk_mutex);
			//	printf("Var x[%d][%d][%d] = %e \n ", subprob, curr_iter, j, x[subprob][curr_iter][j]);
		}
		else if( my_data->sub_data.globals->x[subprob][curr_iter][j] >=
				(1.0-TOLERANCE) ) {
			//	printf("Var x[%d][%d][%d] = %e \n ", subprob, curr_iter, j, x[subprob][curr_iter][j]);
			local_col_name = glp_get_col_name(my_data->sub_data.lp, j);
			i = glp_find_col(original_master_lp, local_col_name );
			my_data->sub_data.globals->final_x[i] =
				/*x[subprob][curr_iter][j]*/ 1.0;
			my_data->sub_data.globals->relaxed_x[i] =
				/*x[subprob][curr_iter][j]*/ 1.0;
			pthread_mutex_lock(&glpk_mutex);
			pre_round_ones++;
			pthread_mutex_unlock(&glpk_mutex);
		}
		else { /* Should never get here. */
			printf("How'd we get here?  x[][][] = %e\n",
					my_data->sub_data.globals->x[subprob][curr_iter][j]);
		}
	}

	pthread_mutex_lock(&glpk_mutex);
	fflush(stdout);
	print_zeros_simple(my_data);
	pthread_mutex_unlock(&glpk_mutex);

	pthread_exit(NULL);
}

