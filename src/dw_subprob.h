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

#ifndef SUBPROB_THREAD_H_
#define SUBPROB_THREAD_H_

#define COMMAND_NULL 0
#define COMMAND_GO   1
#define COMMAND_STOP 2
#define COMMAND_WAIT 3

typedef struct {
	char* infile_name;
	int my_id;
	glp_prob* lp;
	glp_smcp *simplex_control_params;

	//glp_prob* master_read_only;

	int* col_translate;
	/* Other stuff relating to a subproblem... F, c, etc. */
	pthread_mutex_t translate_ready_mutex;
	pthread_mutex_t objective_mutex;
	int translate_ready;

	int local_iteration;

	int phase_one;

	int*    ind;
	double* val;
	int     len;

	double* double_vector;

	int command;

	double** D;
	double*  c;

	double* D_vals;
	int*    D_row_coords;
	int*    D_col_coords;
	int     D_nnz;

	double* condensed_x;

	double* current_solution;

	int num_cols;
	int num_cols_plus;
	double r;
	double obj;
	double new_obj_coeff;

	int unbounded;

	faux_globals* globals;
	master_data* md;
} subprob_struct;

void* subproblem_thread(void* arg);
int signal_availability(subprob_struct* my_data);
void prepare_column(double*, double*, subprob_struct* data);
void organize_solution(subprob_struct*, double*, int);
void printD();

#endif /*SUBPROB_THREAD_H_*/
