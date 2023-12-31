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

#include <semaphore.h>
#include <fcntl.h>
#include <pthread.h>
#include <glpk.h>
#include <stdlib.h>
#include <stdio.h>

#ifndef DECOMPOSE_H_
#define DECOMPOSE_H_

#define BUFF_SIZE             200

#define TOLERANCE             0.000001

#define STACKMEM              2500000
#define DW_INFINITY           1000000.0
#define MAX_PHASE2_ITERATIONS 3000
#define MAX_PHASE1_ITERATIONS 100

/* The 'verbosity' of the output will be determined by the user at launch
 * time.  Default value is OUTPUT_NORMAL.
 */
#define OUTPUT_SILENT   -1
#define OUTPUT_QUIET 	0
#define OUTPUT_NORMAL 	5
#define OUTPUT_VERBOSE 	10
#define OUTPUT_ALL		15

/* Each print statement will include an 'importance' value.  When compared
 * with the 'verbosity' the statement may or may not be executed.  Specifically,
 * if importance <= verbosity, then the message will be printed.
 */
#define IMPORTANCE_DIAG  15
#define IMPORTANCE_MIN   10
#define IMPORTANCE_AVG   5
#define IMPORTANCE_VITAL 0

#define DEFAULT_MIP_GAP  0.01

/* All of these used to be global variables.  I shoved all of them into this
 * struct, created an instance of this struct in main then passed it around
 * one way or another to the functions that needed it.  There are some
 * unfortunate pointer chains as a result in some places, but that is the
 * trade-off, I suppose.
 */
typedef struct {
	double   shift;
	double   mip_gap;
	int      num_clients;
	int      integerize_flag;
	int      rounding_flag;
	int      enforce_sub_integrality;
	int      verbosity;
	int      print_timing_data;
	int      print_final_master;
	int      print_relaxed_sol;
	int      perturb;
	int      max_phase1_iterations;
	int      max_phase2_iterations;
	int      head_service_queue;
	int      tail_service_queue;
	int      write_bases;
	int      write_intermediate_opt_files;
	int      get_monolithic_file;
	int*     service_queue;
	char**   subproblem_names;
	char*    master_name;
	char*    monolithic_name;
	double*  final_x;
	double*  relaxed_x;
	double*** x;
} faux_globals;

/* Globally available variables.  Look into minimizing or eliminating them? */

pthread_attr_t   attr;
pthread_mutex_t  master_lp_ready_mutex;
pthread_cond_t   master_lp_ready_cv;
pthread_mutex_t  service_queue_mutex;
pthread_mutex_t  next_iteration_mutex;
pthread_cond_t   next_iteration_cv;
pthread_mutex_t  master_mutex;
pthread_mutex_t  reduced_cost_mutex;
pthread_mutex_t  glpk_mutex;
pthread_mutex_t  fputs_mutex;
pthread_mutex_t* sub_data_mutex;
#ifdef USE_NAMED_SEMAPHORES
sem_t* customers;
#define CUST_NAMED_SEMAPHORE "/tmp/dw_customers"
#else
sem_t customers;
#endif
glp_prob* original_master_lp;
glp_prob* master_lp;
glp_iocp* parm;
glp_smcp* simplex_control_params;

typedef struct {
	int cols;
	int rows;
	double* row_duals;
	double* c;
} master_data;

typedef struct D_matrix {
	double* values;
	int*    columns;
	int*    pointerB;
	int*    pointerE;
	int     rows;
	int     cols;
	int     rows_plus;
} D_matrix;
D_matrix* D;

typedef struct hook_struct{
	char* outfile_name;
	FILE* outfile;

} hook_struct;

typedef struct new_column{
	char*   name;
	int     non_zeros;
	int*    ind;
	double* val;
} new_column;

typedef struct signal_data{
	int master_lp_ready;
	int next_iteration_ready;
	int current_iteration;
} signal_data;
signal_data* signals;


#endif /*DECOMPOSE_H_*/
