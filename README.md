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

`cmake  ../; make`

Version for debugging is compiled when `-DCMAKE_BUILD_TYPE=Debug` option is
used with cmake; `../` represent the path to the repository (i.e., one
directory up from `build`). The binaries will be in `bin` subdirectory of
`build`. To compile individual, just run `make UTILITY_NAM`

Reference manual is included in the repository (refman.pdf).
