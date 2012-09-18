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
 * dw_blas.c
 *
 *  Created on: Sep 2, 2009
 *      Author: jrios
 *
 *  I wrote these to obviate the need for further libraries.  Granted, you may
 *  get better performance if you use a tested/proven BLAS/SparseBLAS library,
 *  but making the whole project more portable/usable seemed more important.
 *
 *	Search through the code for calls to these functions and you'll see how
 *	simple it would be to use a standard BLAS library call if you wanted to
 *	make the changes.
 */

/* y <- alpha*x + y */
void dw_daxpy(const int len, const double alpha, const double* x, double* y) {
	int i;

	/* Rely on compiler for loop unrolling. */
	for( i = 0; i < len; i++ ) {
		y[i] = alpha*x[i] + y[i];
	}
}

/* return x dot y
 * Note there is no "increment" value for the vectors as there would be in a
 * compliant BLAS implementation.
 */
double dw_ddot(const int len, const double* x, const double* y) {
	double ret = 0.0;
	int i;

	/* Rely on compiler for loop unrolling. */
	for( i = 0; i < len; i++ ) {
		ret += x[i]*y[i];
	}

	return ret;
}

/* return x dot y
 *  where x is in a compressed format.
 * Note there is not increment value for y as there would be in a compliant
 * BLAS implementation.
 */
double dw_ddoti(const int nz, const double* x, const int* ind_x,
		const double* y) {
	double ret = 0.0;
	int i;

	/* Rely on compiler for loop unrolling. */
	for( i = 0; i < nz; i++ ) {
		ret += x[i]*y[ind_x[i]];
	}

	return ret;
}

/* y <- x */
void dw_dcopy(const int len, const double* x, double* y) {
	int i;

	/* Rely on compiler for loop unrolling. */
	for( i = 0; i < len; i++ ) {
		y[i] = x[i];
	}
}

/* A*b = c   where A is stored in coordinate format. */
void dw_dcoogemv(const int m, const int k, const double *val,
		const int *indx, const int *jndx, const int nnz, const double *b,
		double *c) {

	int i,j;
	double *temp_c = c;

	/* Zero-out the product vector. */
	for ( i = 0; i < m; i++ ) temp_c[i] = 0;

	/* Indices are arbitrarily ordered.  Use components of product vector as
	 * accumulators.
	 */
	for( j = 0; j < nnz; j++ ) {
		c[indx[j]] +=  b[jndx[j]] * val[j];
	}
}
