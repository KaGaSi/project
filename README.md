# AnalysisTools

A set of utilities to analyse trajectories from particle-based
molecular simulations. The programs work with [vsf/vcf
files](https://github.com/olenz/vtfplugin/wiki/VTF-format) coordinate
files.

Installation
============

AnalysisTools requires `C` and `Fortran` compilers and `cmake`.  The
compilation should be done in a separate directory, such as `<root>/build`
(where `<root>`} is the AnalysisTools root directory). The following
generates a `Makefile` within `build`:

`mkdir <root>/build; cd <root>/build; cmake ../`

Then run `make` to compile all utilities or `make <utility name>` to
compile a single utility.

The binaries are located in the `<root>build/bin` directory.

Reference manual is included in the repository (refman.pdf).
