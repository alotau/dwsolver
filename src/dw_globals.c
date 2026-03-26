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
 *   SPDX-License-Identifier: GPL-3.0-or-later
 *
 **************************************************************************** */

/*
 * dw_globals.c — single translation unit that *defines* all global variables
 * declared as extern in dw.h.  Every other .c file that includes dw.h sees
 * only the extern declarations and therefore does not produce a second
 * definition.  This resolves the multiple-definition linker errors (KD-001)
 * that prevented the project from building with GCC on Linux.
 *
 * C99 §6.9 ¶5: "an external definition is a declaration that is also a
 * definition of a function or object".  Only one translation unit may provide
 * the external definition for any given identifier with external linkage.
 */

#include "dw.h"

/* Thread synchronisation primitives */
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

/* Semaphore — named on macOS, unnamed on Linux */
#ifdef USE_NAMED_SEMAPHORES
sem_t* customers;
#else
sem_t  customers;
#endif

/* GLPK LP objects and parameters */
glp_prob* original_master_lp;
glp_prob* master_lp;
glp_iocp* parm;
glp_smcp* simplex_control_params;

/* Constraint matrix and signal state */
D_matrix*    D;
signal_data* signals;
