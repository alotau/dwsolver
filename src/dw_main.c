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
 * dantzig_wolfe.c
 *
 * This program implements the Dantzig-Wolfe decomposition algorithm using the
 * GNU Linear Programming Kit.  It is multi-threaded to take advantage of
 * multi-core systems.
 *
 * The most useful reference for this implementation was:
 *
 * Bertsimas and Tsitsiklis. Introduction to Linear Optimization. 1997.
 *
 * Other useful references:
 *
 * Dantzig. Linear Programming and Extensions. 1963.
 * Dantzig and Wolfe. Decomposition Principle for Linear Programs. Operations
 * 		Research. 8(1):101-111. Jan-Feb. 1960.
 * Downey. The Little Book of Semaphores. 2008. (For synchronization)
 * Lasdon. Optimization Theory for Large Systems. 1970.
 *
 * First journal publication using this code:
 *
 * Rios and Ross, Massively Parallel Dantzig-Wolfe Decomposition Applied to
 * 		Traffic Flow Scheduling, Journal of Aerospace Computing, Information,
 * 		and Communication. 7(1):32-45. Jan 2010.
 *
 * Issues:
 * * Relies on hacked version of GLPK wherein the reentrant problems of GLPK
 *    are "fixed."
 * * malloc'ing should be revamped.  Not used in a consistent manner.
 * * sprintf usage is replaced with snprintf in some places, but not all.  And
 *    in places where snprintf is used, the return code isn't always checked.
 * * There are probably several stray variables that aren't useful or are
 *    redundant.
 * * Subproblems do one extra iteration before actually killing selves?
 * * There are many globals that probably don't need to be globals.
 * * DW algorithm is tuned for my problem and isn't tested for more general
 *    inputs.  Some comments are included in the code to point this out.
 * * DW algorithm implemented here doesn't deal with 'infinite' subproblem
 *    solutions.  Master assumes all columns that are offered are bounded
 *    solutions to the subproblem.  This is likely an easy fix within
 *    phase_2_iteration().
 *
 * author: Joseph Rios, Joseph.L.Rios@nasa.gov
 * started: 3/26/08
 *
 *
 */

#include <stdio.h>       /* printf, etc */
#include <stdlib.h>	     /* malloc, etc */
#include <string.h>	     /* strlen, strcmp, etc. */
#include <pthread.h>	 /* Threading */
#include <time.h>        /* For program timing. */
#include "dw.h"
#include "dw_subprob.h"
#include "dw_support.h"
#include "dw_rounding.h"
#include "dw_phases.h"

static int hook(void* info, const char* s);

/* See process_cmdline() in support_functions.c for parameters passed to main.
 */
int main(int argc, char* argv[]) {

	int rc, i, j, num_rows, count, len, num_clients;
	int need_phase_one = 0;
	int row_type;
	int num_phase1_vars;
	int* ind;
	int* col_deletion_indicies = NULL;

	double* val;
	double  variable_value;
	double  perturbation;

	char* local_buffer = (char*) malloc(sizeof(char)*BUFF_SIZE);
	char** phase1_vars;

	void *status;
	FILE* opt_outfile;
	pthread_t* threads;
	subprob_struct* sub_data;

	master_data*  md      = malloc(sizeof(master_data));
	faux_globals* globals = malloc(sizeof(faux_globals));

	/* For clocking the program run time. */
	time_t  t0 = time(NULL);
	clock_t c0 = clock();

	/* Initialize some globals. */
	init_globals(globals);

	/* Process the command line. */
	rc = process_cmdline(argc, argv, globals);
	if( rc == 1 ) {
		dw_printf(IMPORTANCE_VITAL,
				"Something wrong with command line, exiting.\n");
		return 1;
	}
	else if( rc == 2 ) {
		return 0;
	}
	else if( rc != 0 ) {
		return rc;
	}

	dw_printf(IMPORTANCE_AVG,"@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n");
	dw_printf(IMPORTANCE_AVG,"@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n");
	dw_printf(IMPORTANCE_AVG,"DWSOLVER: Stand-alone Dantzig-Wolfe Decomposition Solver\n");
	dw_printf(IMPORTANCE_AVG,"(C) 2010 National Aeronautics and Space Administration\n");
	dw_printf(IMPORTANCE_AVG,"Covered under GPLv3 with Additional Terms\n");
	dw_printf(IMPORTANCE_AVG,"Compiled using GLPK version %s (slightly modified)\n", glp_version());
	dw_printf(IMPORTANCE_AVG,"@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n");
	dw_printf(IMPORTANCE_AVG,"@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n");
	fflush(stdout);

	num_clients   = globals->num_clients; /* Doesn't change. Make local copy.*/
	globals->x    = malloc(sizeof(double**)*num_clients);
	threads       = malloc(sizeof(pthread_t)*num_clients);
	sub_data      = malloc(sizeof(subprob_struct)*num_clients);

	/* Initialize synchronization variables.  */
	init_signals(globals);

	/* Initialize global pthread data structures. */
	init_pthread_data(globals);

	/* Open the output file for the optimization routines */
	strcpy(local_buffer, "out_terminal");
	if( (opt_outfile = fopen(local_buffer, "w")) == NULL) {
		fprintf(stderr, "PROBLEM OPENING TERMINAL OUTFILE. BAILING.\n");
		exit(1);
	}

	/* Force terminal output for optimization routines to file instead */
	glp_term_hook(hook, opt_outfile);

/*****************
 *
 * LAUNCH SUBPROBLEM THREADS
 *
 *****************/

	/* Create/launch the subproblem threads. */
	for( i = 0; i < num_clients; i++) {
		//pthread_mutex_init(&sub_data_mutex[i], NULL); /* done in init_pthread_data() */
		sub_data[i].infile_name = (char*) malloc(sizeof(char)*strlen(globals->subproblem_names[i])+1);
		//printf("%d: Name = %s has len = %d\n", i, globals->subproblem_names[i], strlen(globals->subproblem_names[i]));
		strcpy(sub_data[i].infile_name, globals->subproblem_names[i]);

		/* Initialize subproblem data. */
		sub_data[i].my_id            = i;
		sub_data[i].command          = COMMAND_WAIT;
		sub_data[i].local_iteration  = 0;
		sub_data[i].r                = DW_INFINITY;
		sub_data[i].obj              = DW_INFINITY;
		sub_data[i].phase_one        = 0;
		sub_data[i].globals          = globals;
		sub_data[i].md               = md;

		globals->x[i]           = malloc(sizeof(double*)*MAX_PHASE2_ITERATIONS);

		/* Actually create the thread now. */
		rc = pthread_create(&threads[i],
				&attr,
				subproblem_thread,
				(void *)&sub_data[i]);
		if (rc) {
			dw_printf(IMPORTANCE_VITAL,
					"ERROR; return code from pthread_create() is %d\n", rc);
			exit(-1);
		}
	}

/*****************
 *
 * SET UP MASTER
 *
 *****************/
	/* Open the file containing the original master problem. */
	pthread_mutex_lock(&glpk_mutex);
	original_master_lp = lpx_read_cpxlp(globals->master_name);
	pthread_mutex_unlock(&glpk_mutex);

	/* Prepare the parameters for the simplex method. */
	glp_init_smcp(simplex_control_params);
	simplex_control_params->msg_lev  = GLP_MSG_ON;
	/* Presolving doesn't make sense because many constraints aren't even a
	 * part of the master problem.  Many important vars could be presolved
	 * away. */
	simplex_control_params->presolve = GLP_OFF;

	/* For integer optimization. */
	glp_init_iocp(parm);

	num_rows    = glp_get_num_rows(original_master_lp);

	fflush(stdout);
	/* These arrays will be re-used over and over in GLPK calls. */
	val   = (double*) malloc(sizeof(double)*
			(glp_get_num_cols(original_master_lp) + 1));
	ind   = (int*)    malloc(sizeof(int)*
			(glp_get_num_cols(original_master_lp) + 1));

	/* The 'D' is the part of the constraint matrix for the master that has
	 * to do with the original problem's constraints.  All threads access
	 * this data.  Get it ready.
	 */
	prepare_D(num_rows, ind, val);

	/* Prepare the master data structure used by all threads. */
	prepare_md(md);

	/* Make it easier to find rows and columns by name. */
	glp_create_index(original_master_lp);

	/* Signal all subproblems that master is set up and ready. */
	pthread_mutex_lock(&master_lp_ready_mutex);
	signals->master_lp_ready = 1;
	pthread_cond_broadcast(&master_lp_ready_cv);
	pthread_mutex_unlock(&master_lp_ready_mutex);

	/* Create the _real_ master problem while subprobs do their initial solve.*/
	master_lp = glp_create_prob();
	glp_set_prob_name(master_lp, "master");
	glp_set_obj_dir(master_lp, GLP_MIN);

/*****************
 *
 * ADD ROWS TO MASTER
 *
 *****************/

	/* We have 1 row for each constraint in original master + 1 row for each
	 * subproblem in the form of convexity constraints. */
	glp_add_rows(master_lp, D->rows + num_clients);

	/* Give each convexity row a name and and fix its bounds at 1.0. */
	for( i = (D->rows + 1); i <= (D->rows + num_clients); i++ ) {
		sprintf(local_buffer, "sub%d_convexity", i - D->rows);
		glp_set_row_name(master_lp, i, local_buffer);

		// Fixing row bounds to zero as a first step toward handling bounding
		// rays returned by the subproblems.  Currently dwsolver doesn't
		// handle unbounded subproblems, but if it did, there is a chance
		// there would be no feasible corner points provided by a subproblem
		// and only bounding rays.  If that's the case, the convexity row
		// for the subproblem would have no variables.  If the bound was "1"
		// (as it should be to indicate convexity), then this constraint
		// would be unsatisfiable.  As soon as a 'lambda' gets added for a
		// subproblem, this bound should be set properly to '1'.
		glp_set_row_bnds(master_lp, i, GLP_FX, 1.0, 1.0);
	}

	/* Set the rows relating to original problem.  Give a name based on the
	 * original problem and fix their bounds at b value from original prob. */
	for( i = 1; i <= D->rows; i++ ) {
		glp_set_row_name(master_lp, i, glp_get_row_name(original_master_lp, i));
		glp_set_row_bnds(master_lp, i,
				glp_get_row_type(original_master_lp, i),
				glp_get_row_lb(original_master_lp, i),
				glp_get_row_ub(original_master_lp, i));
	}

	/* We will now add a 'y'/'s' column for each original row. */
	count  = glp_add_cols(master_lp, D->rows);
	dw_printf(IMPORTANCE_DIAG,"Added %d cols to master_lp\n", D->rows);

	/* Each y or s column has a length of 1 if added to formulation. */
	len = 1;

	/* Create enough space to accommodate an auxiliary variable for each row. */
	/* Would be more efficient to determine the amount of memory dynamically. */
	phase1_vars = malloc(sizeof(char*)*D->rows);

	/* Keep track of how many auxiliary variables are added. */
	num_phase1_vars = 0;

	/* Step through each original row and convert it to standard form:
	 *   Ax = b
	 * Note: 'count' should be valued at 0 at this point.
	 */
	dw_printf(IMPORTANCE_DIAG,"There are %d rows.\n", D->rows);
	for( i = count; i < (D->rows + count); i++ ) {
		/* In this for-loop, any variable that is added will only be part of
		 * the row indicated by the index, i, since it is an
		 * auxiliary/slack/surplus variable.
		 */
		ind[1] = i;

		/* Based on row type, add appropriate variable(s). */
		row_type = glp_get_row_type(master_lp, i);

		/* If this is a fixed row, we need an auxiliary variable and, therefore,
		 * will be required to run a phase I procedure.
		 */
		if( row_type == GLP_FX ) {
			phase1_vars[num_phase1_vars] = malloc(sizeof(char)*BUFF_SIZE);
			rc = snprintf(phase1_vars[num_phase1_vars], BUFF_SIZE-1, "y_%d", i);
			if( rc >= BUFF_SIZE - 1) /* Bail if buffer overflow. */
				buffer_overflow( "y_(row_number)", BUFF_SIZE-1);
			/* Sign of auxiliary variable dependent on bound value. */
			val[1] = glp_get_row_ub(master_lp, i) < 0.0 ? -1.0 : 1.0;
			/* Make it part of objective.  To be driven to zero. */
			glp_set_obj_coef(master_lp, i, 1.0);
			glp_set_col_name(master_lp, i, phase1_vars[num_phase1_vars]);
			glp_set_col_bnds(master_lp, i, GLP_LO, 0.0, 0.0);
			glp_set_mat_col (master_lp, i, len, ind, val);
			need_phase_one = 1;
			num_phase1_vars++;
		}
		/* Upper-bounded row requires a surplus variable. */
		else if( row_type == GLP_UP ){
			/* Convert to standard form by making row fixed ('=', not '<='). */
			glp_set_row_bnds(master_lp, i, GLP_FX,
					glp_get_row_ub(master_lp, i), 0.0);

			/* Create a surplus variable for this row. */
			rc = snprintf(local_buffer, BUFF_SIZE-1, "su_%d", i);
			if( rc >= BUFF_SIZE - 1)
				buffer_overflow( "su_(row_number)", BUFF_SIZE-1);
			val[1] = 1.0;
			glp_set_col_name(master_lp, i, local_buffer);
			glp_set_col_bnds(master_lp, i, GLP_LO, 0.0, 0.0);
			glp_set_mat_col (master_lp, i, len, ind, val);
			glp_set_obj_coef(master_lp, i, 0.0);

			/* Since this row is upper-bounded, a negative bound implies we
			 * need a Phase I procedure to find our initial basis.  Here we
			 * add an artificial variable to be driven out.
			 */
			//if( glp_get_row_ub(master_lp, i) < 0 ) {
				int new_col = glp_add_cols(master_lp, 1);
				phase1_vars[num_phase1_vars] = malloc(sizeof(char)*BUFF_SIZE);
				rc = snprintf(phase1_vars[num_phase1_vars],
						BUFF_SIZE-1, "y_%d", i);
				if( rc >= BUFF_SIZE - 1) /* Bail if buffer overflow. */
					buffer_overflow( "y_(row_number)", BUFF_SIZE-1);
				val[1] = -1.0;
				glp_set_obj_coef(master_lp, new_col, 1.0);
				glp_set_col_name(master_lp, new_col,
						phase1_vars[num_phase1_vars]);
				glp_set_col_bnds(master_lp, new_col, GLP_LO, 0.0, 0.0);
				glp_set_mat_col (master_lp, new_col, len, ind, val);

				need_phase_one = 1;
				num_phase1_vars++;
			//}
		}
		/* Lower-bounded row requires a slack variable.
		 * This branch untested/unfinished.  Needs to mirror the previous
		 * branch.  Maybe missing appropriate sign somewhere.  Need examples
		 * to test.
		 */
		else if( row_type == GLP_LO ) { /* Add slack variable. */

			/* Convert to standard form by making row fixed ('=', not '>='). */
			glp_set_row_bnds(master_lp, i, GLP_FX,
					glp_get_row_lb(master_lp, i), 0.0);

			/* Create a slack variable for this row. */
			rc = snprintf(local_buffer, BUFF_SIZE-1, "sl_%d", i);
			if( rc >= BUFF_SIZE - 1)
				buffer_overflow( "sl_(row_number)", BUFF_SIZE-1);
			val[1] = -1.0; /* This is where the sign is in question. */
			glp_set_col_name(master_lp, i, local_buffer);
			glp_set_col_bnds(master_lp, i, GLP_LO, 0.0, 0.0);
			glp_set_mat_col (master_lp, i, len, ind, val);
			glp_set_obj_coef(master_lp, i, 0.0);

			/* Since this row is lower-bounded, a positive bound implies we
			 * need a Phase I procedure to find our initial basis.  Here we
			 * add an artificial variable to be driven out.
			 */
			//if( glp_get_row_lb(master_lp, i) > 0 ) {
				int new_col = glp_add_cols(master_lp, 1);
				phase1_vars[num_phase1_vars] = malloc(sizeof(char)*BUFF_SIZE);
				rc = snprintf(phase1_vars[num_phase1_vars],
						BUFF_SIZE-1, "y_%d", i);
				if( rc >= BUFF_SIZE - 1) /* Bail if buffer overflow. */
					buffer_overflow( "y_(row_number)", BUFF_SIZE-1);
				val[1] = 1.0;
				glp_set_obj_coef(master_lp, new_col, 1.0);
				glp_set_col_name(master_lp, new_col,
						phase1_vars[num_phase1_vars]);
				glp_set_col_bnds(master_lp, new_col, GLP_LO, 0.0, 0.0);
				glp_set_mat_col (master_lp, new_col, len, ind, val);

				need_phase_one = 1;
				num_phase1_vars++;
			//}
		}
		else {
			dw_printf(IMPORTANCE_VITAL,
					"Row isn't fixed or upper/lower bounded?  That is bad.\n");
			dw_printf(IMPORTANCE_VITAL,
					"I'd expect some sort of crash or nonsense solution.\n");
		}

		if( globals->perturb ) { /* Perturb the bound to avoid degeneracy. */
			perturbation = 0.00000001 * i;
			glp_set_row_bnds(master_lp, i, GLP_FX,
				glp_get_row_ub(master_lp, i) + perturbation, 0.0);
		}
	}
/*****************
 *
 * DONE ADDING ROWS TO MASTER
 *
 *****************/

	/* Set the shift, if any. */
	glp_set_obj_coef(master_lp, 0, globals->shift);

	/* Print current problem, just slack/auxiliary variables.  For debugging.*/

	if( glp_get_num_cols(master_lp) > 0 ) {
		pthread_mutex_lock(&glpk_mutex);
		lpx_write_cpxlp(master_lp, "pre_master.cpxlp");
		pthread_mutex_unlock(&glpk_mutex);
	}

	dw_printf(IMPORTANCE_AVG,
		"The master currently has %d rows and %d columns.\n",
		glp_get_num_rows(master_lp), glp_get_num_cols(master_lp));

/*****************
 *
 * BEGIN PHASE ONE (IF NECESSARY)
 *
 *****************/
	/* If we added auxiliary variables, attempt to drive them out. */
	if( need_phase_one ) {
		/* Set the shift. */
		glp_set_obj_coef(master_lp, 0, 0);

		/* Let clients know that we're in Phase I. */
		for( j = 0; j < num_clients; j++ ) {
			pthread_mutex_lock(&sub_data_mutex[j]);
			sub_data[j].phase_one = 1;
			pthread_mutex_unlock(&sub_data_mutex[j]);
		}
		dw_printf(IMPORTANCE_AVG,
				"### Commencing Phase I for reduced master problem...\n");

		/* These data structures will keep track of the objective terms that
		 * the subproblems provide.  Recall that the auxiliary master problem
		 * we are currently solving only has auxiliary variables in the
		 * objective.  When we drive all of them out, the columns that have
		 * been introduced need their appropriate objective terms.
		 */
		char**  obj_names =
			malloc(sizeof(char*)*num_clients*globals->max_phase1_iterations);
		double* obj_coefs =
			malloc(sizeof(double)*num_clients*globals->max_phase1_iterations);
		int*    obj_count =
			malloc(sizeof(int));
		*obj_count = 0;

		if( globals->write_intermediate_opt_files ) {
			sprintf(local_buffer, "phase1_step_0.cpxlp");
			lpx_write_cpxlp(master_lp, local_buffer);
		}

		/* Actually perform the Dantzig-Wolfe algorithm now. */
		for( j = 0; j < globals->max_phase1_iterations; j++ ) {
			dw_printf(IMPORTANCE_AVG,
					"\n###  Iteration %d of phase I ###\n", j+1);
			dw_printf(IMPORTANCE_DIAG, "master_lp has %d cols.\n",
				glp_get_num_cols(master_lp));


			rc = phase_1_iteration(sub_data, globals, (j ? 0 : 1),
					obj_names, obj_coefs, obj_count, md);

			if( globals->write_intermediate_opt_files ) {
				sprintf(local_buffer, "phase1_step_%d.cpxlp", j+1);
				lpx_write_cpxlp(master_lp, local_buffer);
			}

			if( -TOLERANCE < glp_get_obj_val(master_lp) &&
					glp_get_obj_val(master_lp) < TOLERANCE ) {
				dw_printf(IMPORTANCE_AVG,
						"\n@@@@@@@ Minimized auxiliary problem. @@@@@@@\n");
				break;
			}

			if( rc == 0 ) {
				dw_printf(IMPORTANCE_VITAL,"No new columns added and did not reach zero.\n");
				dw_printf(IMPORTANCE_VITAL,"Problem is infeasible.\n");
				dw_printf(IMPORTANCE_VITAL,"Breaking from phase one.\n");
				break;
			}
		}

		dw_printf(IMPORTANCE_AVG,
				"### Preparing for Phase II for reduced master problem...\n");

		/* Let subproblems know we are out of phase one. */
		for( j = 0; j < num_clients; j++ ) {
			pthread_mutex_lock(&sub_data_mutex[j]);
			sub_data[j].phase_one = 0;
			pthread_mutex_unlock(&sub_data_mutex[j]);
		}

		/* Delete the y columns. */
		col_deletion_indicies = malloc(sizeof(int)*(D->rows + 1));
		for( i = 0; i < num_phase1_vars; i++ ) {
			col_deletion_indicies[i+1] =
				glp_find_col(master_lp, phase1_vars[i]);
			free(phase1_vars[i]);
		}
		glp_del_cols(master_lp, num_phase1_vars, col_deletion_indicies);
		free( col_deletion_indicies ) ;

		/* Restore objective coeff's of the current variables in system. */
		if( globals->verbosity >= OUTPUT_NORMAL )
			printf("Trying to restore %d objective coeff's...\n", *obj_count);
		for( i = 0; i < *obj_count; i++ ) {
			if( glp_find_col(master_lp, obj_names[i]) <= 0 ) {
				printf("Didn't find a column I expected to find.\n");
				printf("Things are probably broken.\n");
			}
			glp_set_obj_coef(master_lp, glp_find_col(master_lp, obj_names[i]),
					obj_coefs[i]);
			free(obj_names[i]);
		}

		/* Rebuild the basis after removing columns. This is necessary in the
		 * presence of degeneracy, but safe to do regardless.
		 */
		glp_adv_basis(master_lp, 0);
		rc = glp_simplex(master_lp, simplex_control_params);
		dw_printf(IMPORTANCE_AVG,
				"Recomputed basis and solved.  Simplex returned %d\n", rc);

		/* Reset the appropriate shift. */
		dw_printf(IMPORTANCE_DIAG, "THE SHIFT IS %3.2f\n", globals->shift);
		glp_set_obj_coef(master_lp, 0, globals->shift);

		/* Clean up. */
		free(obj_count);
		free(obj_coefs);
		free(phase1_vars);
		free(obj_names);
		dw_printf(IMPORTANCE_AVG,"Phase II will start now.\n");
	}
	else {
		dw_printf(IMPORTANCE_AVG,
				"No auxiliary variables introduced.  Straight to Phase II.\n");
		free(phase1_vars);
	}
/*****************
 *
 * END PHASE I
 *
 *****************/

	dw_printf(IMPORTANCE_DIAG,"There will be a maximum of %d DW iterations.\n",
			MAX_PHASE2_ITERATIONS);

/*****************
 *
 * BEGIN PHASE II
 *
 *****************/
	/* Actually perform Dantzig-Wolfe algorithm now. */
	for( j = 0; j < globals->max_phase2_iterations; j++ ) {

		dw_printf(IMPORTANCE_AVG,"\n###  Iteration %d of phase II ###\n",
				j);
		dw_printf(IMPORTANCE_DIAG,"master_lp has %d cols.\n",
				glp_get_num_cols(master_lp));

		if( globals->write_intermediate_opt_files ) {
			sprintf(local_buffer, "master_step_%d.cpxlp", j);
			lpx_write_cpxlp(master_lp, local_buffer);
		}

		rc = phase_2_iteration(sub_data, globals, md);

		/* These lines may be useful to someone solving a relaxed integer
		 * instance.  Or maybe not.
		 */
		//check_degeneracy();
		//if( check_col_integrality() )
		//	printf("Non-integral variables!\n");
		//else printf("All primal values are integral.\n");

		if( rc == 0 ) {
			dw_printf(IMPORTANCE_AVG,"Didn't add any columns?  Let's break.\n");
			for(i = 0; i < num_clients; i++) {
				/* Unfortunately, this command doesn't reach the subproblems
				 * until they've completed their current iteration. On
				 * later iterations, this may not be so bad, but if there is
				 * no guarantee on how long an iteration may take. This can
				 * be easily(?) refactored to send the correct code
				 * dependent upon the value of rc in phase_2_iteration().
				 */
				pthread_mutex_lock(&sub_data_mutex[i]);
				sub_data[i].command = COMMAND_STOP;
				pthread_mutex_unlock(&sub_data_mutex[i]);
			}
			pthread_mutex_lock(&next_iteration_mutex);
			signals->current_iteration++;
			pthread_cond_broadcast(&next_iteration_cv);
			pthread_mutex_unlock(&next_iteration_mutex);
			break;
		}
		else {
			dw_printf(IMPORTANCE_DIAG,"### Master added %d columns.\n", rc);
		}
	}
/*****************
 *
 * END PHASE II
 *
 *****************/

	if( globals->perturb ) {
		for( i = count; i < (D->rows + count); i++ ) {
			perturbation = 0.00000001 * i;
			glp_set_row_bnds(master_lp, i, GLP_FX,
				glp_get_row_ub(master_lp, i) - perturbation, 0.0);
		}
		glp_simplex(master_lp, simplex_control_params);
		if( check_col_integrality() ) {

			dw_printf(IMPORTANCE_AVG,
					"There are basic variables with non-integer values.\n");
			//printf("But I'm not doing anything about it.\n");
			//		printf("I want to see if there is an equivalent integer solution.  Solving...\n");
			//		simplex_rc = glp_intopt(master_lp, parm);
			//		printf("Integer optimization returned %d.\n", simplex_rc);
			//		printf("Integer optimal value is %3.2f.\n", glp_mip_obj_val(master_lp));
		}
		else dw_printf(IMPORTANCE_AVG,"The basic variables are all integer.\n");

		dw_printf(IMPORTANCE_AVG,"#########################################\n");
		dw_printf(IMPORTANCE_AVG,"####  Master objective value = %e \n",
				glp_get_obj_val(master_lp));
		dw_printf(IMPORTANCE_AVG,"#########################################\n");

	}

	/* Print timing before trying integer optimization. */
	dw_printf(IMPORTANCE_AVG,"Done with solving the relaxation...\n");
	if( globals->print_timing_data ) {
		print_timing(t0, c0 );
		t0 = time(NULL);
		c0 = clock();
	}

	/* Print out the final version of the master problem. Solving this on its
	 * own later will provide the optimal value.
	 */
	if( globals->print_final_master ) {
		dw_printf(IMPORTANCE_AVG,"%-40s ","Printing final master to done.cpxlp...");
		glp_write_lp(master_lp, NULL, "done.cpxlp");
		dw_printf(IMPORTANCE_AVG,"DONE!\n");
	}

	/* Wait for subthreads to join. */
	dw_printf(IMPORTANCE_AVG, "%-40s ", "Waiting for subthreads...");
	for(i=0; i<num_clients; i++) {
		rc = pthread_join(threads[i], &status);
		if (rc) {
			dw_printf(IMPORTANCE_VITAL,
					"ERROR; thread %d return code from pthread_join() is %d\n",
					i, rc);
			dw_printf(IMPORTANCE_VITAL,
					"Going to continue waiting for all joins,");
			dw_printf(IMPORTANCE_VITAL,
					" but behavior now unpredictable.\n");
			//exit(-1);
		}
			dw_printf(IMPORTANCE_DIAG,
				"### Completed join with thread %d status = %ld\n",i,
				(long)status);
	}
	dw_printf(IMPORTANCE_AVG,"DONE!\n");

	/* Get the optimal values for each variable from original problem. */
	//get_solution(sub_data);

	/* Print out the final values of the convexity variables.  To make this
	 * more useful, you may want to modify this to dump to a file. */
	if( globals->verbosity >= OUTPUT_ALL )
		for( i = num_rows; i <= glp_get_num_cols(master_lp); i++ ) {
			printf("  %s: \t%3.2f\t(%3.2f)\n", glp_get_col_name(master_lp, i),
					glp_get_col_prim(master_lp, i),
					glp_mip_col_val(master_lp, i));
		}

	/* Print relaxed solution to file... */
	if( globals->print_relaxed_sol ) {
		dw_printf(IMPORTANCE_AVG, "%-40s ","Printing relaxed solution to file...");
		process_solution(sub_data, RELAXED_SOLN_FILE, DW_MODE_PRINT_RELAXED);
		dw_printf(IMPORTANCE_AVG, "DONE!\n");
	}

/*****************
 *
 * MAKE RELAXED SOLUTION INTEGRAL?
 *
 *****************/
	/* Try rounding... */
	if( globals->rounding_flag ) {
		dw_printf(IMPORTANCE_AVG, "%-40s ", "Calculating rounded solution...");
		process_solution(sub_data, ROUNDED_ZEROS_FILE, DW_MODE_ROUND_SOL);
		if( globals->print_timing_data ) {
			print_timing(t0, c0 );
			t0 = time(NULL);
			c0 = clock();
		}
		dw_printf(IMPORTANCE_AVG, "DONE!\n");
	}

	/* Try 'choosing one column per subproblem'... */
	if( globals->integerize_flag ) {
		/* Get an integer solution. */
		dw_printf(IMPORTANCE_AVG, "Going to integerize by enforcing binary constraint on lambdas.\n");
		dw_printf(IMPORTANCE_AVG, "Note that this isn't guaranteed to work well or at all.\n");
		glp_iocp* int_parm = malloc(sizeof(glp_iocp));
		glp_init_iocp(int_parm);
		int_parm->mip_gap = globals->mip_gap;
		glp_intopt(master_lp, int_parm);
		//lpx_intopt(master_lp);
		dw_printf(IMPORTANCE_AVG, "Integer optimal solution: %3.1f\n",
				glp_mip_obj_val(master_lp));
		dw_printf(IMPORTANCE_AVG,
				"Replacing LP relaxed solution with integer solution...\n");
		for( i = 1; i <= glp_get_num_cols(master_lp); i++ ) {
			variable_value = glp_mip_col_val(master_lp, i);
			glp_set_col_bnds(master_lp, i, GLP_FX, variable_value, 0.0);
		}
		dw_printf(IMPORTANCE_AVG, "Setting basis...\n");
		glp_simplex(master_lp, simplex_control_params);

		FILE* integerization_zero_file;
		if( (integerization_zero_file = fopen(INTEGERIZED_ZEROS_FILE, "w")) == NULL ) {
			printf("Problem opening file: %s\n", INTEGERIZED_ZEROS_FILE);
			return -1;
		}

		process_solution(sub_data, INTEGERIZED_ZEROS_FILE, DW_MODE_PRINT_RELAXED);
		//process_solution(sub_data, INTEGERIZED_ZEROS_FILE, DW_MODE_INT_SOL);
//		for( i = 1; i <= glp_get_num_cols(master_lp); i++ ) {
//			parse_zero_var(glp_get_col_prim(master_lp, i), i,
//					master_lp, integerization_zero_file);
//		}

		/* Collect runtime information and print it out. */
		dw_printf(IMPORTANCE_AVG,
				"Done with solving the integer optimization...\n");
		if( globals->print_timing_data ) {
			print_timing(t0, c0 );
			t0 = time(NULL);
			c0 = clock();
		}
	}

/*****************
 *
 * CLEAN UP
 *
 *****************/
	/* Clean up after ourselves. */
	free_sub_data(sub_data, globals);
	free_globals(globals, md);
	free(globals);
	free(ind);
	free(val);
	free(threads);
	free(simplex_control_params);
	free(local_buffer);

	dw_printf(IMPORTANCE_AVG,
			"Master made it to the end. Exiting gracefully, dignity intact.\n");
	return 0;
}

/* The function that glp_term_hook uses to redirect terminal output */
static int hook(void* info, const char* s) {
	FILE* outfile = info;
	pthread_mutex_lock(&fputs_mutex);
	fputs(s, outfile);
	pthread_mutex_unlock(&fputs_mutex);
	return 1;
}

