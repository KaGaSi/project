# AnalysisTools

A bundle of linux programs to analyse trajectories from coarse-grained
molecular dynamics simulations. The programs work with
[vsf/vcf files](https://github.com/olenz/vtfplugin/wiki/VTF-format)
coordinate files.

Installation
===

All programs can be compiled using cmake, which generates Makefile, and
subsequently make. It requires C and FORTRAN compilers (I use gcc and
gcc-fortran). It's better to compile the package in its own directory.

For example, create directory `build`; `cd` to the directory; and
to install all the programs, run from the command line this
command to first create unix makefile and then compile all programs:

`cmake -DCMAKE_BUILD_TYPE=Release -G "Unix Makefiles" ../ ; make`

Version for debugging is compiled when `-DCMAKE_BUILD_TYPE=Debug` is used
instead of `-DCMAKE_BUILD_TYPE=Release` and `../` represent the path to the
repository (i.e., one directory up from `build`). The binaries will be in
`bin` subdirectory of `build`.

To compile individual C programs using `gcc`, run:

`gcc -O3 -lm $PATH_TO_REPOSITORY/AnalysisTools.c $PATH_TO_REPOSITORY/Options.c $PATH_TO_PROGRAM -o $OUTPUT_NAME`

Reference manual is included in the repository (refman.pdf), but the
documentation can be generated anew simply by running `make` in the
repository directory. Aside from refman.pdf, it also creates `doc/html`
directory that contains html documentation. All documentation was created
using [doxygen](http://www.stack.nl/~dimitri/doxygen/).
