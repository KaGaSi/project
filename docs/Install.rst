Installation
============

AnalysisTools require ``C`` and ``Fortran`` compilers as well as ``cmake``.
It is recommended to build the utilities in a separate directory (e.g.,
``build``). Assuming the current working directory is the AnalysTools root
directory, the following generates a ``Makefile`` within ``build``::

   mkdir build; cd build; cmake ../

Then to compile all utilities, just run ``make`` in the ``build``
directory. Individual utilities can be compiled by running ``make
<utility>`` instead.

The binaries are located in the ``build/bin`` directory.
