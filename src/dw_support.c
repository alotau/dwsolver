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
 * support_functions.c
 *
 *  Created on: Jun 3, 2009
 *      Author: jrios
 */
#include "dw.h"
#include "dw_subprob.h"
#include "dw_support.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdarg.h>
#include <pthread.h>
#include <config.h>
#include <glpk.h>

int dw_verbosity = OUTPUT_NORMAL;

void dw_printf(int importance_level, char* format_string, ...) {
	/* See dw.h for documentation on the interaction between verbosity
	 * and importance.
	 */
	va_list vl;

	if( importance_level <= dw_verbosity ) {
		va_start(vl, format_string);
		vprintf(format_string, vl);
		va_end(vl);
	}
}

/* This just looks at all of the basic variables and sees which, if any, have
 * a 0.0 primal value.  If any basic variables are valued 0.0, this implies
 * degeneracy.  No big deal.
 */
void check_degeneracy() {
	int i;
	for( i = 1; i <= glp_get_num_cols(master_lp); i++ ) {
		if( glp_get_col_stat(master_lp, i) == GLP_BS ) {
			if( glp_get_col_prim(master_lp, i) < TOLERANCE ) {
				dw_printf(IMPORTANCE_DIAG,
						"Column %d (%s) is basic and has value %3.2f\n",
						i, glp_get_col_name(master_lp, i),
						glp_get_col_prim(master_lp, i));
			}
		}
		else if( glp_get_col_stat(master_lp, i) == GLP_NU ) {
			dw_printf(IMPORTANCE_DIAG,
					"Column %d (%s) is GLP_NU and has value %3.2f\n",
					i, glp_get_col_name(master_lp, i),
					glp_get_col_prim(master_lp, i));
		}
	}
}

/* The 'D' matrix is essentially the coefficeint matrix for the coupling
 * constraints.  These constraints are read from the 'master' input file.
 * This function is called by the master thread and prepares the D matrix
 * for easy global access by the subthreads.  The subthreads use this matrix
 * to offer columns to the master problem. */
void prepare_D(int num_rows, int* ind, double* val) {
	int count, i, j, len;
	double sign;

	D           = malloc(sizeof(D_matrix));
	D->cols     = glp_get_num_cols(original_master_lp);
	D->rows     = num_rows;
	D->rows_plus= num_rows+1;
	D->values   =
		(double*) malloc(sizeof(double)*glp_get_num_nz(original_master_lp));
	D->columns  = (int*) malloc(sizeof(int)*glp_get_num_nz(original_master_lp));
	D->pointerB = (int*) malloc(sizeof(int)*num_rows);
	D->pointerE = (int*) malloc(sizeof(int)*num_rows);

	count = 0;  /* 1 for 1-based indexing, 0 for 0-based indexing? */
	for( i = 1; i <= num_rows; i++ ) {
		sign = 1.0;
		if( glp_get_row_type(original_master_lp, i) == GLP_LO ) {
			if (glp_get_row_lb(original_master_lp, i) < 0.0 ) sign = -1.0;
		}
		else if( glp_get_row_type(original_master_lp, i) == GLP_UP ) {
			if( glp_get_row_ub(original_master_lp, i) < 0.0 ) sign = -1.0;
		}
		D->pointerB[i-1] = count;
		len = glp_get_mat_row(original_master_lp, i, ind, val);
		for( j = 1; j <= len; j++ ) {
			D->values[count]  = sign * val[j];
			D->columns[count] = ind[j]-1;
			count++;
		}
		D->pointerE[i-1] = count;
	}
}

/* Prepare master data.  Most of this is static data based on the 'master'
 * input file.  row_duals is updated each iteration with the dual cost
 * of each coupling constraint.  This is used by the subproblems to set
 * their objective functions.  */
void prepare_md(master_data* md) {
	int i;
	md->row_duals = (double*) calloc((D->rows+1), sizeof(double));
	md->c         = (double*) calloc((D->cols+1), sizeof(double));
	/* Store the objective coefficients for the original problem. */
	for( i = 1; i <= D->cols; i++ ) {
		md->c[i] = glp_get_obj_coef(original_master_lp, i);
	}
	md->cols = glp_get_num_cols(original_master_lp);
	md->rows = glp_get_num_rows(original_master_lp);
}

/* Initialize global variables.  Do all of these NEED to be in here?  Hard to
 * say anymore.  Better than being true globals?  Hard to say as well. */
void init_globals(faux_globals* fg) {
	fg->verbosity = OUTPUT_NORMAL;
	fg->integerize_flag = 0;
	fg->perturb = 0;
	fg->get_monolithic_file = 1;
	fg->max_phase1_iterations = MAX_PHASE1_ITERATIONS;
	fg->max_phase2_iterations = MAX_PHASE2_ITERATIONS;
	fg->shift = 0.0;
	fg->print_timing_data = 0;
	fg->print_final_master = 1;
	fg->print_relaxed_sol = 1;
	fg->write_bases = 0;
	fg->write_intermediate_opt_files = 0;
	fg->rounding_flag = 0;
	fg->enforce_sub_integrality = 0;
	simplex_control_params = (glp_smcp*) malloc(sizeof(glp_smcp));
	parm = malloc(sizeof(glp_iocp));

	/* Global variables for all threads to put/pop to/from service queue. */
	fg->head_service_queue = 0;
	fg->tail_service_queue = 0;
}

/* These are basically flags that get set/reset throughout algorithm. */
void init_signals(faux_globals* fg) {
	fg->service_queue = (int*) malloc(sizeof(int)*fg->num_clients);
	signals       = (signal_data*) malloc(sizeof(signal_data));

	/* Initialize synchronization variables.  */
	signals->current_iteration    = 0;
	signals->master_lp_ready      = 0;
	signals->next_iteration_ready = 0;
}

/* Initialize global pthread data structures. */
void init_pthread_data(faux_globals* fg) {
	size_t stacksize;
	int rc, i;
	pthread_mutexattr_t* my_mutex_attr = malloc(sizeof(pthread_mutexattr_t));

	/* Initialize mutual exclusion and condition variable objects */
	pthread_mutexattr_init(my_mutex_attr);

	sub_data_mutex =
		(pthread_mutex_t*) malloc(sizeof(pthread_mutex_t)*fg->num_clients);
	for( i = 0; i < fg->num_clients; i++ )
		pthread_mutex_init(&(sub_data_mutex[i]), NULL);
	pthread_mutex_init(&master_lp_ready_mutex, NULL);
	pthread_cond_init (&master_lp_ready_cv, NULL);
	pthread_mutex_init(&service_queue_mutex, NULL);
	pthread_mutex_init(&next_iteration_mutex, my_mutex_attr);
	pthread_mutex_init(&master_mutex, NULL);
	pthread_mutex_init(&reduced_cost_mutex, NULL);
	pthread_mutex_init(&glpk_mutex, NULL);
	pthread_mutex_init(&fputs_mutex, NULL);
	pthread_cond_init (&next_iteration_cv, NULL);
	pthread_attr_init (&attr);
#ifdef USE_NAMED_SEMAPHORES
	customers = sem_open(CUST_NAMED_SEMAPHORE, O_CREAT, S_IRWXU, 0);
#else
	sem_init(&customers, 0, 0);
#endif
	/* Make sure each thread is join-able. */
	pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

	/* Give each thread plenty of stack memory.  Probably unnecessary. */
	pthread_attr_getstacksize (&attr, &stacksize);
	stacksize = sizeof(double)*STACKMEM;
	rc = pthread_attr_setstacksize (&attr, stacksize);
	pthread_attr_getstacksize (&attr, &stacksize);

	free(my_mutex_attr);
}

/* This is just a clean-up function.  Trying to be tidy with memory. */
void free_globals(faux_globals* fg, master_data* md) {

	int i;
	dw_printf(IMPORTANCE_AVG, "%-40s ","Freeing globals...");
	/* Free D matrix */
	free(D->columns);
	free(D->values);
	free(D->pointerB);
	free(D->pointerE);
	free(D);

	/* Free main glp objects. */
	glp_delete_prob(original_master_lp);
	glp_delete_prob(master_lp);


	/* Free file names. */
	for( i = 0; i < fg->num_clients; i++ ) {
		free(fg->subproblem_names[i]);
		pthread_mutex_destroy(&sub_data_mutex[i]);
	}
	free(fg->subproblem_names);
	free(fg->master_name);
	free(fg->monolithic_name);

	/* Free misc. stuff. */
	free(fg->service_queue);
	free(md->row_duals);
	free(md->c);
	free(md);
	free(sub_data_mutex);
	free(fg->x);
	free(signals);

	/* Clean up mutual exclusion stuff. */
	pthread_mutex_destroy(&service_queue_mutex);
	pthread_mutex_destroy(&master_lp_ready_mutex);
#ifdef USE_NAMED_SEMAPHORES
	sem_close(customers);
	sem_unlink(CUST_NAMED_SEMAPHORE);
#else
	sem_destroy(&customers);
#endif

	pthread_attr_destroy(&attr);

	dw_printf(IMPORTANCE_AVG, "DONE!\n");
}

/* Purging non-basic variables should work well in theory.  In practice, I've
 * found that it increases runtime and can lead to some cycling.  Maybe there
 * is a better way to implement this function? Not currently used anywhere.
 */
void purge_nonbasics() {
	int i;
	int count = 0;
	int* nb_indices = malloc(sizeof(int)* (glp_get_num_cols(master_lp)+1));

	dw_printf(IMPORTANCE_AVG, "There are currently %d structural variables.\n",
			glp_get_num_cols(master_lp));
	dw_printf(IMPORTANCE_AVG, "There are currently %d rows.\n",
			glp_get_num_rows(master_lp));

	/* Check each column.  If non-basic, schedule for deletion by storing
	 * index in nb_indices.  Spare the surplus variables from deletion.
	 */
	for( i = 1 ; i <= glp_get_num_cols(master_lp); i++ ) {
		if( glp_get_col_name(master_lp, i)[0] == 's' ) ;
		else if( glp_get_col_stat(master_lp, i) != GLP_BS ) {
			count++;
			nb_indices[count] = i;
		}
	}

	dw_printf(IMPORTANCE_AVG, "Purging %d non-basic variables...\n", count);

	glp_del_cols(master_lp, count, nb_indices);

	glp_adv_basis(master_lp, 0);

	free(nb_indices);
}

void print_usage(int argc, char* argv[]) {
	int i;
	printf("Usage: %s -g <guide_file_name> [options]\n\n", argv[0]);
	printf("This program was called with the following command line:\n");
	for( i = 0; i < argc; i++ ) {
		printf("%s ", argv[i]);
	}
	printf("\n\n");
	printf(" -v, --version           Print version information and exit.\n");
	printf(" -g, --guidefile <file>  REQUIRED. Use <file> as the guidefile.\n");
	printf(" -h, -help, --help       Print this usage and return.\n");
	printf(" -i, --integerize        After solving the LP with Dantzig-Wolfe, try\n");
	printf("                         some integerization of the solution by\n");
	printf("                         enforcing integrality constraints on the \n");
	printf("                         convexity variables. This may be computationally\n");
	printf("                         intensive.  This destroys relaxed solution.\n");
	printf(" -e, --sub-int-enforce   Enforce integer solutions in the subproblems at\n");
	printf("                         each iteration. If using -i you may want/need this\n");
	printf("                         flag as well.\n");
	printf(" --mip-gap <f>           Set the mip gap to f. Only used if -i used as well.\n", DEFAULT_MIP_GAP);
	printf("                         Default %f.\n", DEFAULT_MIP_GAP);
	printf(" -r, --round             After solving the LP with Dantzig-Wofle, round\n");
	printf("                         all convexity variables to nearest int. This will\n");
	printf("                         likely break some constraints, but is a fast method\n");
	printf("                         for obtaining an integer \"solution\".\n");
	printf(" -p, --perturb           Perturb RHS (Experimental).\n");
	printf(" --verbose               Verbose mode.\n");
	printf(" --output-all            For diagnostic purposes only. More output than verbose mode.\n");
	printf(" -n, --output-normal     Normal output mode (default)\n");
	printf(" -q, --quiet             Quiet mode. (Still some minimal output.)\n");
	printf(" --silent                Silent mode. (Maybe not completely silent?)\n");
	printf(" --skip-monolithic-read  The guide file does not contain a monolithic file.\n");
	printf(" --write-bases           Write out the basis information at each iteration.\n");
	printf(" --write-int-probs       Write out intermediate optimization files.\n");
	printf(" --write-final-master    Write out the final master problem (default).\n");
	printf(" --no-write-final-master Do not write out the final master problem.\n");
	printf(" --print-timing          Print runtime information. \n");
	printf("                         Note: Ignores print level flags, i.e. will try to print no matter what.\n");
	printf(" --no-timing             Do not write out timing information (default).\n");
	printf(" --phase1_max <n>        Allow at most n phase I iterations. Default %d.\n", MAX_PHASE1_ITERATIONS);
	printf(" --phase2_max <n>        Allow at most n phase II iterations. Default %d.\n", MAX_PHASE2_ITERATIONS);
	printf("\n");
	printf("At a minimum, this program requires a guidefile.\n");
	printf("The guidefile will direct the program to read the appropriate\n");
	printf("data files.  The guide file will have the following format:\n\n");
	printf(" n\n");
	printf(" file_1\n");
	printf(" file_2\n");
	printf(" ...\n");
	printf(" file_n\n");
	printf(" master_file\n");
	printf(" monolithic_file\n");
	printf(" objective_constant\n");
	printf(" <eof>\n\n");
	printf("All of the listed files are in CPLEX format.  n is the number\n");
	printf("of subproblems.  monolithic_file is essentially a file with all\n");
	printf("of the constraints/variables contained in all of the subproblem\n");
	printf("files and the master_file.  Variable names must be consistent\n");
	printf("throughout all files.  The objective_constant is optional.\n");
	printf("\n");
}

/* Simple processing of command line. Probably many more options would be
 * useful. */
int process_cmdline(int argc, char* argv[], faux_globals* fg) {
 	int i;
 	int opened_guide_file = 0;
 	FILE* input_file = NULL;
 	FILE** subproblem_files;
 	FILE* master_file;
 	FILE* monolithic_file;
 	char* buffer = (char*) malloc(sizeof(char)*BUFF_SIZE);

 	for( i = 1; i < argc; i++ ) {
 		/* Specifying the guide file. */
 		if( strcmp(argv[i], "-g") == 0 || strcmp(argv[i], "--guidefile") == 0) {
 			i++;
 			if( (input_file = fopen(argv[i], "r")) == NULL ) {
 				dw_printf(IMPORTANCE_VITAL,
 							"Problem opening the guide file: %s\n", argv[i]);
 				print_usage(argc, argv);
 				return 1;
 			}
 			else opened_guide_file = 1;
 		}
 		/* Help message. */
 		else if( strcmp(argv[i], "-h") == 0 || strcmp(argv[i], "--help") == 0 ||
 				strcmp(argv[i], "-help") == 0) {
 			print_usage(argc, argv);
 			return 2;
 		}
 		else if( strcmp(argv[i], "-v") == 0 ||
 				strcmp(argv[i], "--version") == 0) {
 			dw_printf(IMPORTANCE_AVG, "%s: Dantzig-Wolfe Solver version %s\n",
 					PACKAGE, VERSION);
 			dw_printf(IMPORTANCE_AVG, "Built using a modified GNU Linear Programming Kit, version %s\n\n", glp_version());
 			dw_printf(IMPORTANCE_AVG, "Copyright (C) 2010 Joseph Rios, National Aeronautics and Space Administration.\n");
 			dw_printf(IMPORTANCE_AVG, "All rights reserved.\n\n");
 			dw_printf(IMPORTANCE_AVG, "This program has absolutely no warranty.\n\n");
 			dw_printf(IMPORTANCE_AVG, "This program is free software; you may re-distribute it under the terms\n");
 			dw_printf(IMPORTANCE_AVG, "of the GNU General Public License version 3 or later.\n");
 			return 2;
 		}
 		else if(strcmp(argv[i], "--output-all") == 0) {
 			fg->verbosity = OUTPUT_ALL;
 		}
 		else if(strcmp(argv[i], "--verbose") == 0) {
 			fg->verbosity = OUTPUT_VERBOSE;
 		}
 		else if( strcmp(argv[i], "-q") == 0 ||
 				strcmp(argv[i], "--quiet") == 0) {
 			fg->verbosity = OUTPUT_QUIET;
 		}
 		else if( strcmp(argv[i], "-n") == 0 ||
 				strcmp(argv[i], "--normal-output") == 0) {
 			fg->verbosity = OUTPUT_NORMAL;
 		}
 		else if( strcmp(argv[i], "-p") == 0 ||
 				strcmp(argv[i], "--perturb") == 0) {
 			fg->perturb = 1;
 		}
 		else if( strcmp(argv[i], "--print-timing") == 0 ) {
 			fg->print_timing_data = 1;
 		}
 		else if( strcmp(argv[i], "--no-timing") == 0 ) {
 			fg->print_timing_data = 0;
 		}
 		else if( strcmp(argv[i], "--write-final-master") == 0 ) {
 			fg->print_final_master = 1;
 		}
 		else if( strcmp(argv[i], "--no-write-final-master") == 0 ) {
 			fg->print_final_master = 0;
 		}
 		else if( strcmp(argv[i], "-i") == 0 ||
 				strcmp(argv[i], "--integerize") == 0 ) {
 			fg->integerize_flag = 1;
 		}
 		else if( strcmp(argv[i], "-e") == 0 ||
 				strcmp(argv[i], "sub-int-enforce") == 0) {
 			fg->enforce_sub_integrality = 1;
 		}
 		else if( strcmp(argv[i], "--mip_gap") == 0 ||
 				 strcmp(argv[i], "--mip-gap") == 0) {
 			if( ((fg->mip_gap = atof (argv[++i]))) <= 0 )	{
 				dw_printf(IMPORTANCE_VITAL,
 						"%s is not a valid number for mip_gap. Ignoring.\n", argv[i]);
 				fg->mip_gap = DEFAULT_MIP_GAP;
 			}
 		}
 		else if( strcmp(argv[i], "-r") == 0
 				|| strcmp(argv[i], "--round") == 0 ) {
 			fg->rounding_flag = 1;
 		}
 		else if( strcmp(argv[i], "--write-bases") == 0 ) {
 			fg->write_bases = 1;
 		}
 		else if( strcmp(argv[i], "--write-int-probs") == 0 ) {
 			fg->write_intermediate_opt_files = 1;
 		}
 		else if( strcmp(argv[i], "--skip-monolithic-read") == 0 ) {
 			fg->get_monolithic_file = 0;
 		}
 		else if( strcmp(argv[i], "--phase1_max") == 0 ) {
 		 	if( ((fg->max_phase1_iterations = atoi (argv[++i]))) <= 0 )	{
 		 		dw_printf(IMPORTANCE_VITAL,
 		 					"%s is not a valid number of iterations. Ignoring.\n", argv[i]);
 		 		fg->max_phase1_iterations = MAX_PHASE1_ITERATIONS;
 		 	}
 		}
 		else if( strcmp(argv[i], "--phase2_max") == 0 ) {
 			if( ((fg->max_phase2_iterations = atoi (argv[++i]))) <= 0 )	{
 				dw_printf(IMPORTANCE_VITAL,
 							"%s is not a valid number of iterations. Ignoring.\n", argv[i]);
 				fg->max_phase2_iterations = MAX_PHASE2_ITERATIONS;
 			}
 		}
 		else printf("Command line argument %s not recognized. Trying to continue.\n", argv[i]);
 	}
 	if( !opened_guide_file ) {
 		print_usage(argc, argv);
 		return 1;
 	}

 	fgets(buffer, BUFF_SIZE, input_file);
 	fg->num_clients = atoi(buffer);
 	if( fg->num_clients < 1 || fg->num_clients > 26000 ) {
 		dw_printf(IMPORTANCE_VITAL,
 				"Exiting: There is an strange number of subproblems: %d\n",
 				fg->num_clients);
 		return 1;
 	}
 	else {
 		dw_printf(IMPORTANCE_DIAG,
 				"I think there are %d subproblems.\n", fg->num_clients);
 	}

 	fg->subproblem_names = malloc(sizeof(char*)*fg->num_clients);
 	subproblem_files     = malloc(sizeof(FILE*)*fg->num_clients);
 	fg->master_name      = malloc(sizeof(char)*BUFF_SIZE);
 	fg->monolithic_name  = malloc(sizeof(char)*BUFF_SIZE);

 	/* Get number of clients, then get each client file's name. */
 	for( i = 0; i < fg->num_clients; i++ ) {
 		fg->subproblem_names[i] = malloc(sizeof(char)*BUFF_SIZE);
 		fgets(fg->subproblem_names[i], BUFF_SIZE, input_file);
 		fg->subproblem_names[i][strlen(fg->subproblem_names[i])-1] = '\0';
 		if( (subproblem_files[i] = fopen(fg->subproblem_names[i], "r")) == NULL ) {
 			dw_printf(IMPORTANCE_VITAL, "Problem opening file: %s\n", fg->subproblem_names[i]);
 			return 1;
 		}
 		else {
 			if( fg->verbosity >= OUTPUT_VERBOSE ) {
 				dw_printf(IMPORTANCE_DIAG,
 						"Opened %s for reading.\n", fg->subproblem_names[i]);
 				dw_printf(IMPORTANCE_DIAG, "Now closing it...\n");
 			}
 			fclose(subproblem_files[i]);
 		}
 	}

 	/* Get master file name. */
 	fgets(fg->master_name, BUFF_SIZE, input_file);
 	fg->master_name[strlen(fg->master_name)-1] = '\0';
 	if( (master_file = fopen(fg->master_name, "r")) == NULL ) {
 		dw_printf(IMPORTANCE_VITAL, "Problem opening master file: %s\n", fg->master_name);
 		return 1;
 	}
 	else {
 		dw_printf(IMPORTANCE_DIAG, "Opened %s for reading.\n", fg->master_name);
 		dw_printf(IMPORTANCE_DIAG, "Now closing it.\n");
 		fclose(master_file);
 	}

 	/* Get monolithic file's name. */
 	/* This code is vestigial.  Leaving it in for backward compatibility with
 	 * my previously generated guidefiles. Should remove this feature
 	 * altogether and fix my old guidefiles...
 	 */
 	if( fg->get_monolithic_file ) {
 		fgets(fg->monolithic_name, BUFF_SIZE, input_file);
 		fg->monolithic_name[strlen(fg->monolithic_name)-1] = '\0';
 		/*
 		if( (monolithic_file = fopen(fg->monolithic_name, "r")) == NULL ) {
 			printf("Problem opening master file: %s\n", fg->monolithic_name);
 			return 1;
 		}
 		else {

 			if( fg->verbosity >= OUTPUT_VERBOSE ) {
 				printf("Opened %s for reading.\n", fg->monolithic_name);
 				printf("Now closing it.\n");
 			}

 			fclose(monolithic_file);

 		}
 		*/
 	}
 	else fg->monolithic_name = '\0';

 	/* Get the "shift" if it's included in the guide file. */
 	strcpy(buffer, "\0"); /* "Clear" the buffer. */
 	fgets(buffer, BUFF_SIZE, input_file);
 	if( buffer != NULL ) {
 		fg->shift = atof(buffer);
 	}
 	/* else fg->shift stays set to 0.0 as initialized. */

 	dw_verbosity = fg->verbosity;

 	free(buffer);
 	free(subproblem_files);
 	return 0;
}

/* This is a VERY problem-specific function that no one is really going to be
 * able to use.  It parses a very specifically-formatted variable name.
 * Should be generalized somehow or excluded?
 */
int parse_zero_var(double value, int index, glp_prob* lp, FILE* zero_file) {
	const char* var_name;
	char* local_col_name;
	char* sector_name;

	if( value < TOLERANCE) {
		var_name = glp_get_col_name(lp, index);
		local_col_name = malloc(sizeof(char)*BUFF_SIZE);
		sector_name = malloc(sizeof(char)*BUFF_SIZE);
		strcpy(local_col_name, var_name);
		if( strtok(local_col_name, "(") == NULL ) printf("NULL 1");
		strcpy( sector_name, strtok(NULL, ",") );
		fprintf(zero_file, "%s %s\n", sector_name, strtok(NULL, ","));
		free(sector_name);
		free(local_col_name);
	}

	return 0;
}

void write_basis(int iteration) {
	int i;
	int basic_var_count = 0;
	char* filename = malloc(sizeof(char)*BUFF_SIZE);
	FILE* basis_file;
	sprintf(filename, "basis_iteration_%d", iteration);
	if( (basis_file = fopen(filename, "w")) == NULL 	) {
		dw_printf(IMPORTANCE_VITAL, "Problem opening %s for writing.\n",
				filename);
		return;
	}

	for( i = 1; i <= glp_get_num_cols(master_lp); i++ ) {
		if( glp_get_col_stat(master_lp, i) == GLP_BS ) {
			fprintf(basis_file, "%s value %.8e\n", glp_get_col_name(master_lp, i), glp_get_col_prim(master_lp, i));
			basic_var_count++;
		}
	}
	fprintf(basis_file, "Num basics: %d\n", basic_var_count);
	fprintf(basis_file, "Num cols  : %d\n", glp_get_num_cols(master_lp));
	free(filename);
	fclose(basis_file);
}

/* Checks to see if the reduced master's variables are integral (binary) or
 * not.  Useful information if you are applying Dantzig-Wolfe to an integer
 * program and are enforcing integrality in the subproblems.
 */
int check_col_integrality() {
	int i;
	int rc = 0;
	for( i = 1; i <= glp_get_num_cols(master_lp); i++ ) {
		if( glp_get_col_prim(master_lp, i) >= TOLERANCE &&
			glp_get_col_prim(master_lp, i) <= 1.0-TOLERANCE){
			rc++;
		}
	}
	return rc;
}

void check_aux_vars() {
	int i;
	for( i = 1; i <= D->rows; i++ ) {
		if( glp_get_col_prim(master_lp, i) <= TOLERANCE ) {
			/* This is actually would be "IMPORTANCE_DIAG" information, but
			 * we assume the user wouldn't have called this function if he/she
			 * didn't want to see it printed.
			 */
			dw_printf(IMPORTANCE_AVG, "   %s: %s = %e\n", glp_get_row_name(master_lp, i),
					glp_get_col_name(master_lp, i),
					glp_get_col_prim(master_lp, i));
		}
	}
}
/* Print to the screen the values of the variables from the original problem
 * as discovered by the DW algorithm.  Expects a specifcally-formatted variable
 * name.  Should generalize this function.
 */
void get_solution(subprob_struct* sub_data) {
	int i;
	int j;
	int subprob;
	int iteration;
	const char* col_name;
	char* local_col_name = malloc(sizeof(char)*BUFF_SIZE);
	for( i = D->rows_plus; i <= glp_get_num_cols(master_lp); i++ ) {
		if( glp_get_col_prim(master_lp, i) != 0.0 ) {
			col_name = glp_get_col_name(master_lp, i);
			strcpy(local_col_name, col_name);
			/* Tease out the subprob number and the iteration. */
			strtok(local_col_name, "_");
			/* Would probably make sense to do safety checks in here,
			 * but for now just take the chance. */
			subprob = atoi(strtok(NULL, "_"));
			iteration = atoi(strtok(NULL, "_"));
			/* This is actually would be "IMPORTANCE_DIAG" information, but
			 * we assume the user wouldn't have called this function if he/she
			 * didn't want to see it printed.
			 */
			dw_printf(IMPORTANCE_AVG, "%s: %d, %d\n", col_name, subprob, iteration);

			for( j = 1; j < sub_data[subprob].num_cols_plus; j++ ) {
				if( glp_get_col_prim(master_lp, i) != 1.0) {
					dw_printf(IMPORTANCE_AVG, "%s:", glp_get_col_name(sub_data[subprob].lp, j));
					dw_printf(IMPORTANCE_AVG, "    %3.2f\n",
							sub_data->globals->x[subprob][iteration][j]);
				}
			}
		}
	}
	free( local_col_name ) ;
}

/* A quick and dirty check for initial feasibility.  This isn't sufficient
 * enough to declare a problem feasible, but it may find a trivially
 * infeasible constraint.
 */
int dirty_feas_check() {
	int i, j;
	int num_cols;
	int sum;
	int bad = 0;
	double bound;
	double* row_coef = malloc(sizeof(double)*glp_get_num_cols(original_master_lp));
	for( i = 1; i <= glp_get_num_rows(original_master_lp); i++) {
		dw_printf(IMPORTANCE_AVG, "Row %d... ", i);
		sum = 0;
		bound = glp_get_row_ub(original_master_lp, i);
		if( bound >= 0.0 ) {
			dw_printf(IMPORTANCE_AVG, "bound is %.2f, skipping the check.\n", bound);
			continue;
		}
		num_cols = glp_get_mat_row(original_master_lp, i, NULL, row_coef);

		for( j = 1; j <= num_cols; j++) {
			if( bound < 0.0 && row_coef[j] == -1.0 ) sum--;
			//else if( bound > 0.0 && row_coef[j] == 1.0) sum++;
		}

		if( bound < 0.0 && sum > bound ) {
			dw_printf(IMPORTANCE_AVG, "Constraint %d (%s) has upper bound of ",
					i, glp_get_row_name(original_master_lp, i));
			dw_printf(IMPORTANCE_AVG,
					"%.2f but can only get to %d.\n", bound, sum);
			bad++;
		}
		else dw_printf(IMPORTANCE_AVG,
				"this may be feasible... bound = %.2f, sum = %d\n", bound, sum);


	}
	if( bad > 0 ) printf("This problem is infeasible.\n");
	else printf("This problem may be feasible.\n");

	return bad;
}

void buffer_overflow(char* str, int len) {
	dw_printf(IMPORTANCE_VITAL, "Whoa. I just tried to overflow a buffer.  The");
	dw_printf(IMPORTANCE_VITAL, " program is likely unstable now. I will bail.\n");
	dw_printf(IMPORTANCE_VITAL, "Tried to write: %s into buffer of size %d.\n",
			str, len);
	exit(1);
}


/* Collect runtime information and print it out. */
void print_timing( time_t start_time, clock_t start_clock ) {

	printf("$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$\n");
	clock_t c1 = clock();
	time_t  t1 = time(NULL);
	printf ("\n TIMING INFORMATION:\n\n");
	printf ("\t begin  (CPU time):       %d\n", (int) start_clock);
	printf ("\t end    (CPU time):       %d\n", (int) c1);
	printf ("\t elapsed CPU time:        %f s\n",
			(float) (c1 - start_clock)/CLOCKS_PER_SEC);
	printf ("\n");
	printf ("\t begin (wall clock):      %ld\n", (long) start_time);
	printf ("\t end   (wall clock):      %ld\n", (long) t1);
	printf ("\t elapsed wall clock time: %ld\n", (long) (t1 - start_time));

	printf("$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$\n\n");
}

/* Frees all memory associated with the subprob_struct data. */
void free_sub_data(subprob_struct* sub_data, faux_globals* fg) {
	int i, j;

	dw_printf(IMPORTANCE_AVG,"%-40s ",
			"Freeing subproblem data..");
	for(i = 0; i < fg->num_clients; i++ ) {
		dw_printf(IMPORTANCE_DIAG, "\n   Freeing sub[%d] data...  ", i);
		free(sub_data[i].c);
		free(sub_data[i].col_translate);
		free(sub_data[i].ind);
		free(sub_data[i].val);
		free(sub_data[i].double_vector);
		free(sub_data[i].infile_name);
		free(sub_data[i].condensed_x);
		free(sub_data[i].D_col_coords);
		free(sub_data[i].D_row_coords);
		free(sub_data[i].D_vals);
		for( j = 0; j < signals->current_iteration - 1; j++) {
			if( fg->x[i][j] != NULL ) free(fg->x[i][j]);
		}
		free(fg->x[i]);
		glp_delete_prob(sub_data[i].lp);
	}
	free(sub_data);

	dw_printf(IMPORTANCE_AVG, "DONE!\n");
}
