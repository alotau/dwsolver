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
 * dw_rounding.h
 *
 *  Created on: Dec 15, 2008
 *      Author: jrios
 */

#ifndef ROUNDING_H_
#define ROUNDING_H_

#define DW_MODE_PRINT_RELAXED 91
#define DW_MODE_ROUND_SOL     92
#define DW_MODE_INT_SOL       93

typedef struct solution_info sol_info;

#define ROUNDED_ZEROS_FILE     "zeros_rounded"
#define INTEGERIZED_ZEROS_FILE "integerized_zeros"
#define RELAXED_SOLN_FILE      "relaxed_solution"

struct solution_info {
	int subprob;
	int iteration;
	double value;
	sol_info* next;
};

typedef struct {
	sol_info* si;
	int num_si;
	int id;
	FILE* zero_file;
	subprob_struct sub_data;
} int_thread_data;

void* rounding_thread(void* arg) ;
void* solution_printing_thread(void* arg) ;
void check_broken_constraints();
void print_zeros(int_thread_data*);
void print_zeros_simple(int_thread_data*);
int round_solution(subprob_struct* sub_data, char* zero_file_name);
int process_solution(subprob_struct* sub_data, char* out_file_name, int mode);

#endif /* ROUNDING_H_ */
