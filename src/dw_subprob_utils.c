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
#include "dw_subprob_utils.h"

void initialize_subproblem(subprob_struct* my_data, glp_smcp* simplex_control_params) {
    int id = my_data->my_id;
    my_data->simplex_control_params = simplex_control_params;

    // Set up GLPK parameters
    glp_init_smcp(simplex_control_params);
    simplex_control_params->presolve = GLP_OFF;
    simplex_control_params->msg_lev  = GLP_MSG_ERR;

    // Copy global data
    my_data->globals->verbosity = my_data->globals->verbosity;

    dw_printf(IMPORTANCE_DIAG, "Subproblem %d initialized.\n", id);
}

glp_prob* setup_and_solve_initial_problem(subprob_struct* my_data, glp_smcp* simplex_control_params) {
    pthread_mutex_lock(&glpk_mutex);
    glp_prob* lp = lpx_read_cpxlp(my_data->infile_name);
    pthread_mutex_unlock(&glpk_mutex);

    glp_iocp* int_parm = malloc(sizeof(glp_iocp));
    glp_init_iocp(int_parm);
    int_parm->msg_lev = GLP_MSG_ERR;

    for (int i = 1; i <= glp_get_num_cols(lp); i++) {
        if (glp_get_col_type(lp, i) == GLP_LO) {
            glp_set_col_bnds(lp, i, GLP_DB, glp_get_col_lb(lp, i), 3000.0);
        }
    }

    int ret = glp_simplex(lp, simplex_control_params);
    my_data->obj = glp_get_obj_val(lp);

    if (my_data->globals->enforce_sub_integrality) {
        glp_intopt(lp, int_parm);
        my_data->obj = glp_mip_obj_val(lp);
    }

    free(int_parm);
    return lp;
}

void wait_for_master_ready() {
    pthread_mutex_lock(&master_lp_ready_mutex);
    if (!signals->master_lp_ready) {
        pthread_cond_wait(&master_lp_ready_cv, &master_lp_ready_mutex);
    }
    pthread_mutex_unlock(&master_lp_ready_mutex);
}

void setup_data_structures(subprob_struct* my_data, glp_prob* lp) {
    int my_col = glp_get_num_cols(lp) + 1;
    my_data->num_cols = glp_get_num_cols(lp);
    my_data->num_cols_plus = my_col;

    my_data->col_translate = malloc(sizeof(int) * my_col);
    my_data->c = calloc(my_col, sizeof(double));
    my_data->condensed_x = malloc(sizeof(double) * my_col);
    my_data->D_vals = malloc(sizeof(double));
    my_data->D_row_coords = malloc(sizeof(int));
    my_data->D_col_coords = malloc(sizeof(int));
    my_data->D_nnz = 0;

    dw_printf(IMPORTANCE_DIAG, "Data structures set up for subproblem %d.\n", my_data->my_id);
}