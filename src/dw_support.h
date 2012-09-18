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
 * support_functions.h
 *
 *  Created on: Jun 3, 2009
 *      Author: jrios
 */

#ifndef SUPPORT_FUNCTIONS_H_
#define SUPPORT_FUNCTIONS_H_

//struct faux_globals;

void solve(glp_prob *prob, glp_smcp *params, const char* solution_file);
int process_cmdline(int argc, char* argv[], faux_globals*);
void free_globals(faux_globals* fg, master_data* md);
//int phase_1_iteration(subprob_struct* sub_data, int, char**, double*, int*);
//int phase_2_iteration(subprob_struct* sub_data);
void dw_printf(int,char*,...);
void test_matrix_math();
void get_solution();
void print_timing( time_t start_time, clock_t start_clock ) ;
void free_sub_data(subprob_struct* sub_data, faux_globals*) ;
int check_col_integrality() ;
void check_aux_vars();
void write_basis(int);
void check_degeneracy();
void purge_nonbasics();
int dirty_feas_check();
void init_pthread_data(faux_globals*);
void init_globals(faux_globals*);
void init_signals(faux_globals*);
void prepare_D(int, int*, double*);
void prepare_md(master_data* md);
void buffer_overflow(char* str, int len);
int parse_zero_var(double value, int index, glp_prob* lp, FILE* zero_file);

#endif /* SUPPORT_FUNCTIONS_H_ */
