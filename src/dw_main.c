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
 * dw_main.c — Thin CLI wrapper around libdwsolver.
 *
 * Parses command-line arguments and a guide file, then delegates all solver
 * work to dw_solve() from the callable library.
 */

#include "dw_solver.h"
#include <config.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#define CLI_BUFF_SIZE 200

/* Strip trailing \r and \n. */
static void cli_chomp(char *s) {
	size_t len = strlen(s);
	while (len > 0 && (s[len-1] == '\n' || s[len-1] == '\r'))
		s[--len] = '\0';
}

static void print_usage(int argc, char *argv[]) {
	int i;
	printf("Usage: %s -g <guide_file_name> [options]\n\n", argv[0]);
	printf("This program was called with the following command line:\n");
	for( i = 0; i < argc; i++ ) printf("%s ", argv[i]);
	printf("\n\n");
	printf(" -v, --version           Print version information and exit.\n");
	printf(" -g, --guidefile <file>  REQUIRED. Use <file> as the guidefile.\n");
	printf(" -h, -help, --help       Print this usage and return.\n");
	printf(" -i, --integerize        Integerize the solution after LP solve.\n");
	printf(" -e, --sub-int-enforce   Enforce integer subproblem solutions.\n");
	printf(" --mip-gap <f>           Set the MIP gap (default %g).\n", 0.01);
	printf(" -r, --round             Round convexity variables to nearest int.\n");
	printf(" -p, --perturb           Perturb RHS (experimental).\n");
	printf(" --verbose               Verbose output.\n");
	printf(" --output-all            Maximum diagnostic output.\n");
	printf(" -n, --output-normal     Normal output (default).\n");
	printf(" -q, --quiet             Quiet output.\n");
	printf(" --silent                Silent mode.\n");
	printf(" --skip-monolithic-read  Guide file has no monolithic entry.\n");
	printf(" --write-bases           Write bases each iteration.\n");
	printf(" --write-int-probs       Write intermediate LP files.\n");
	printf(" --write-final-master    Write final master LP (default off).\n");
	printf(" --no-write-final-master Do not write final master LP.\n");
	printf(" --print-timing          Print runtime information.\n");
	printf(" --no-timing             Suppress timing output.\n");
	printf(" --phase1_max <n>        Max phase I iterations (default 100).\n");
	printf(" --phase2_max <n>        Max phase II iterations (default 3000).\n");
	printf("\n");
	printf("The guide file format:\n");
	printf("  n\n  sub_file_1\n  ...\n  sub_file_n\n  master_file\n");
	printf("  monolithic_file\n  objective_constant\n\n");
}

/*
 * parse_cli — parse argv into opts and populate file-path out-parameters.
 *
 * Returns 0 on success, 1 on error, 2 on informational exit (--help/-v).
 * On success, *master_out and *subs_out are heap-allocated; caller must free.
 */
static int parse_cli(int argc, char *argv[],
                     dw_options_t *opts,
                     char **master_out,
                     char ***subs_out,
                     int *nsubs_out) {
	int i, n;
	int skip_mono = 0;
	int opened_guide = 0;
	FILE *gf = NULL;
	char buf[CLI_BUFF_SIZE];
	char **subs = NULL;

	for( i = 1; i < argc; i++ ) {
		if( strcmp(argv[i], "-g") == 0 || strcmp(argv[i], "--guidefile") == 0 ) {
			if( i + 1 >= argc ) {
				fprintf(stderr, "Option %s requires a guide file name.\n", argv[i]);
				print_usage(argc, argv);
				return 1;
			}
			i++;
			gf = fopen(argv[i], "r");
			if( !gf ) {
				fprintf(stderr, "Cannot open guide file: %s\n", argv[i]);
				print_usage(argc, argv);
				return 1;
			}
			opened_guide = 1;
		}
		else if( strcmp(argv[i], "-h") == 0 || strcmp(argv[i], "--help") == 0 ||
		         strcmp(argv[i], "-help") == 0 ) {
			print_usage(argc, argv);
			return 2;
		}
		else if( strcmp(argv[i], "-v") == 0 || strcmp(argv[i], "--version") == 0 ) {
			printf("%s\n", dw_version());
			return 2;
		}
		else if( strcmp(argv[i], "--output-all") == 0 )    opts->verbosity = DW_OUTPUT_ALL;
		else if( strcmp(argv[i], "--verbose") == 0 )       opts->verbosity = DW_OUTPUT_VERBOSE;
		else if( strcmp(argv[i], "-n") == 0 ||
		         strcmp(argv[i], "--output-normal") == 0 ) opts->verbosity = DW_OUTPUT_NORMAL;
		else if( strcmp(argv[i], "-q") == 0 ||
		         strcmp(argv[i], "--quiet") == 0 )         opts->verbosity = DW_OUTPUT_QUIET;
		else if( strcmp(argv[i], "--silent") == 0 )        opts->verbosity = DW_OUTPUT_SILENT;
		else if( strcmp(argv[i], "-p") == 0 ||
		         strcmp(argv[i], "--perturb") == 0 )       opts->perturb = 1;
		else if( strcmp(argv[i], "--print-timing") == 0 )  opts->print_timing_data = 1;
		else if( strcmp(argv[i], "--no-timing") == 0 )     opts->print_timing_data = 0;
		else if( strcmp(argv[i], "--write-final-master") == 0 )    opts->print_final_master = 1;
		else if( strcmp(argv[i], "--no-write-final-master") == 0 ) opts->print_final_master = 0;
		else if( strcmp(argv[i], "-i") == 0 ||
		         strcmp(argv[i], "--integerize") == 0 )    opts->integerize_flag = 1;
		else if( strcmp(argv[i], "-e") == 0 ||
		         strcmp(argv[i], "--sub-int-enforce") == 0 ) opts->enforce_sub_integrality = 1;
		else if( strcmp(argv[i], "--mip_gap") == 0 ||
		         strcmp(argv[i], "--mip-gap") == 0 ) {
			if( i + 1 >= argc ) {
				fprintf(stderr, "Option %s requires a value.\n", argv[i]);
				print_usage(argc, argv);
				return 1;
			}
			double v = atof(argv[++i]);
			if( v > 0.0 ) opts->mip_gap = v;
			else fprintf(stderr, "%s is not a valid mip_gap; using default.\n", argv[i]);
		}
		else if( strcmp(argv[i], "-r") == 0 ||
		         strcmp(argv[i], "--round") == 0 )         opts->rounding_flag = 1;
		else if( strcmp(argv[i], "--write-bases") == 0 )
			fprintf(stderr, "Warning: --write-bases is not supported by this build; option ignored.\n");
		else if( strcmp(argv[i], "--write-int-probs") == 0 )
			fprintf(stderr, "Warning: --write-int-probs is not supported by this build; option ignored.\n");
		else if( strcmp(argv[i], "--skip-monolithic-read") == 0 ) skip_mono = 1;
		else if( strcmp(argv[i], "--phase1_max") == 0 ) {
			if( i + 1 >= argc ) {
				fprintf(stderr, "Option %s requires a value.\n", argv[i]);
				print_usage(argc, argv);
				return 1;
			}
			int v = atoi(argv[++i]);
			if( v > 0 ) opts->max_phase1_iterations = v;
			else fprintf(stderr, "%s is not a valid iteration count; using default.\n", argv[i]);
		}
		else if( strcmp(argv[i], "--phase2_max") == 0 ) {
			if( i + 1 >= argc ) {
				fprintf(stderr, "Option %s requires a value.\n", argv[i]);
				print_usage(argc, argv);
				return 1;
			}
			int v = atoi(argv[++i]);
			if( v > 0 ) opts->max_phase2_iterations = v;
			else fprintf(stderr, "%s is not a valid iteration count; using default.\n", argv[i]);
		}
		else printf("Command line argument %s not recognized. Trying to continue.\n", argv[i]);
	}

	if( !opened_guide ) {
		print_usage(argc, argv);
		return 1;
	}

	/* Read number of subproblems. */
	if( !fgets(buf, CLI_BUFF_SIZE, gf) ) {
		fprintf(stderr, "Guide file is empty.\n"); fclose(gf); return 1;
	}
	n = atoi(buf);
	if( n < 1 || n > 26000 ) {
		fprintf(stderr, "Invalid subproblem count: %d\n", n); fclose(gf); return 1;
	}

	subs = malloc(sizeof(char*) * (size_t)n);
	if( !subs ) { fclose(gf); return 1; }

	/* Read subproblem file names. */
	for( i = 0; i < n; i++ ) {
		subs[i] = malloc(CLI_BUFF_SIZE);
		if( !subs[i] ) { fclose(gf); return 1; }
		if( !fgets(subs[i], CLI_BUFF_SIZE, gf) ) {
			fprintf(stderr, "Guide file truncated reading subproblem %d\n", i);
			fclose(gf); return 1;
		}
		cli_chomp(subs[i]);
	}

	/* Read master file name. */
	*master_out = malloc(CLI_BUFF_SIZE);
	if( !*master_out ) { fclose(gf); return 1; }
	if( !fgets(*master_out, CLI_BUFF_SIZE, gf) ) {
		fprintf(stderr, "Guide file truncated reading master file name\n");
		fclose(gf); return 1;
	}
	cli_chomp(*master_out);

	/* Read (and discard) monolithic file name for guide-file backward compat. */
	if( !skip_mono ) {
		if( fgets(buf, CLI_BUFF_SIZE, gf) ) { /* consume line, ignore */ }
	}

	/* Read objective shift (optional). */
	buf[0] = '\0';
	if( fgets(buf, CLI_BUFF_SIZE, gf) && buf[0] != '\0' ) {
		double s = atof(buf);
		if( s != 0.0 ) opts->shift = s;
	}

	fclose(gf);
	*subs_out  = subs;
	*nsubs_out = n;
	return 0;
}

/* -------------------------------------------------------------------------
 * main — thin CLI wrapper
 * ------------------------------------------------------------------------- */
int main(int argc, char *argv[]) {
	dw_options_t  opts;
	dw_result_t   result = {0, 0.0, 0, NULL};
	char         *master_file    = NULL;
	char        **subproblem_files = NULL;
	int           num_subproblems  = 0;
	int           rc, i;

	dw_options_init(&opts);

	rc = parse_cli(argc, argv, &opts, &master_file, &subproblem_files, &num_subproblems);
	if( rc == 1 ) {
		fprintf(stderr, "Something wrong with command line, exiting.\n");
		return 1;
	}
	else if( rc == 2 ) {
		return 0;
	}
	else if( rc != 0 ) {
		return rc;
	}

	rc = dw_solve(master_file,
	              (const char *const *)subproblem_files,
	              num_subproblems,
	              &opts,
	              &result);

	/* Release file-path allocations from parse_cli. */
	for( i = 0; i < num_subproblems; i++ )
		free(subproblem_files[i]);
	free(subproblem_files);
	free(master_file);

	dw_result_free(&result);
	return rc;
}


