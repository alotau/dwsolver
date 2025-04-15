#ifndef DW_SUBPROB_UTILS_H
#define DW_SUBPROB_UTILS_H

#include <glpk.h>
#include "dw.h"

void initialize_subproblem(subprob_struct* my_data, glp_smcp* simplex_control_params);
glp_prob* setup_and_solve_initial_problem(subprob_struct* my_data, glp_smcp* simplex_control_params);
void wait_for_master_ready();
void setup_data_structures(subprob_struct* my_data, glp_prob* lp);
void solve_subproblem_iterations(subprob_struct* my_data, glp_prob* lp, glp_smcp* simplex_control_params);

#endif // DW_SUBPROB_UTILS_H