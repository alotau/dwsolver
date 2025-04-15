# DWSOLVER (Dantzig-Wolfe Solver)

## Overview
The Dantzig-Wolfe Solver program is a stand-alone implementation of the
Dantzig-Wolfe Decomposition algorithm.  The GNU Linear Programming Kit provides
the functions for all of the necessary linear programming (reading in problems,
performing the simplex algorithm, querying various LP data strucures, etc.).
This software also uses POSIX Threads (pthreads) for parallel solving of 
subproblems.  The implementation is general in that any linear program that can
be presented to the software in block-angular form (with some current 
limitations) can be solved.

## Known Limitations
This repo was updated to take advantage of GH Actions over a decade after it was 
previously touched. I could only get the MacOS runner to work properly. Using the
Ubuntu runner gives errors at the linking phase due to some loose usage of variable
definitions within [dw.h](./src/dw.h). Trying a few things to fix old C code like
trying to make vars static, using externs, etc. This was a fools errand. Trying to update multithreaded C code from years ago may be the most dangerous game known to programmers. That said, I would happily take any PRs that show this code being built by a GH Action on an Ubuntu runner.

## Examples and Tests
There is a an [examples](./examples/) directory with several problems plucked from popular textbooks and some websites. There is also a toy version of the problem that initiated the writing of this code related to air traffic management. Each example has a README and should run successfully with a properly built dwsolver executable. There is also a new [tests](./tests/) directory that contains a [test script](./tests/dw-tests.sh) that runs most of the examples.  A couple of the examples have non-deterministic solutions (but deterministic optimum values), and I didn't take the time to write tests that check for the correct optimum for those problems. Again, happy for any PR that runs those examples, parses the output files to pluck out the optimum, and shows it is the expected value.

## Copyright
Copyright  2010 United States Government National Aeronautics and Space 
Administration (NASA).  No copyright is claimed in the United States under 
Title 17, U.S. Code. All Other Rights Reserved.

## Legal
This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

In accordance with the GNU General Public License Version 3, 29 June 2007
(GPL V3) Section 7. Additional Terms, Additional Permissions are added as
exceptions to the terms of the GPL V3 for this program.  These additional
terms should have been received with this program in a file entitled
["ADDITIONAL_LICENSE_TERMS"](./ADDITIONAL_LICENSE_TERMS).  If a copy was not provided, you may request
one from the contact author listed below.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

## Additional Info
See the file [COPYING](COPYING) for the GNU Public License.

See the file [INSTALL](INSTALL) for compilation and installation instructions.

All Dantzig-Wolfe source files created by the author begins with the prefix
"dw" while all other source files in the src/ directory are GLPK files.  All
GLPK files are provided as published by the author(s) of GLPK except for those
specified in the .patch file provided with this release.  The changes to GLPK
were made to implement a thread-friendly version of GLPK necessary for 
implementation of a parallel Dantzig-Wolfe algorithm.

~~Please report bugs/comments/suggestions/patches to Joseph.L.Rios@nasa.gov.~~
I'm not longer with NASA, so best to just submit through GitHub as an issue.
