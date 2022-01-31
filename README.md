# AnalysisTools

A bundle of utilities to analyse trajectories from particle-based molecular
simulations. The programs work mainly with [vsf/vcf
files](https://github.com/olenz/vtfplugin/wiki/VTF-format) trajectory files.

Installation
===

All programs can be compiled using cmake. It requires C and FORTRAN compilers
(I use gcc and gcc-fortran). The following is a simple oneliner to compile the
utilities in its own directory (assuming one is in the AnalysisTools root
directory):

`mkdir build; cd build; cmake ../; make`

The compiled binaries are then in the `build/bin` directory.

Reference manual is included in the repository (manual.pdf).
