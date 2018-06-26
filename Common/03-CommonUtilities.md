# Common utilities {#Common}

Utilies that are not specific to any given system and are used for all
simulations.

# Aggregates & Aggregates-NotSameBeads utility {#Aggregates}

These utilities determine which molecules belong to which aggregates
according to a simple criterion: two molecules belong to the same aggregate
if they share at least a specified number of contact pairs. A contact pair
is a pair of two beads belonging to different molecules which are closer
than certain distance. Both the distance and the number of needed contact
pairs are arguments of the command as well as bead types to consider.
Specified molecule(s) can be excluded from aggregate calculation (both from
aggregate calculation and the output `.agg` file).

While the Aggregates utility uses all possible pairs of given bead types,
Aggregates-NotSameBeads does not use same-type pairs. For example, if bead types
`A` and `B` are given, Aggregates will use all three bead type pairs,
that is `A-A`, `A-B` and `B-B` (provided the beads are not in the same molecules), but
Aggregates-NotSameBeads will not use `A-A` or `B-B`. Therefore, at least two
bead types must be provided in Aggregates-NotSameBeads argument.

Usage:

`Aggregates` (or `Aggregates-NotSameBeads`) ` <input.vcf> <distance> <contacts> <output.agg> <type names> <options>`

> `<input.vcf>`
> > input coordinate filename (must end with `.vcf`) containing indexed timesteps
> `<distance>`
> > minimum distance for two beads to be in contact (constituting one
> > contact pair)
> `<contacts>`
> > minimum number of contact pairs to consider two molecules to be in one
> > aggregate
> `<output.agg>`
> > output filename (must end with `.agg`) containing information about
> > aggregates
> `<type names>`
> > names of bead types to use for calculating contact pairs
> `<options>`
> > `-x <name(s)>`
> > > exclude specified molecule(s) from calculation of aggregates
> > `-j <joined.vcf>`
> > > filename for coordinates of joined aggregates (must end with `.vcf`)

# Average utility {#Average}

Average uses binning method to analyse data stored in a supplied file. It
prints average, statistical error and estimate of integrated
autocorrelation time (tau). Empty lines and lines beginning with '#' are
skipped. The program prints to the screen four numbers: `<n_blocks>
<simple average> <statistical error> <estimate of tau>`.

A way to estimate
a 'real' value of tau is to use a wide range of `<n_blocks>` and then plot
`<tau>` as a function of `<n_blocks>`. Since the number of data points in a
block has to be larger than tau (say, 10 times larger), plotting
`<number of data lines>/10/<n_blocks>` vs. `<n_blocks>` will produce an
exponential function that will intersect the plotted `<tau>`. A value of tau
near the intersection (but to the left where the exponential is above `<tau>`)
can be considered a good estimate for tau.

Usage:

`Average <filename> <column> <discard> <n_blocks>`

> `<filename>`
> > name of data filel
> `<column>`
> > column number in the file containing the data to analyze
> `<discard>`
> > number of data values considered as equilibrium
> `<n_blocks>`
> > number of blocks for binning analysis

# BondLength utility — not extensively used or tested {#BondLength}

BondLength utility calculates normalised distribution of bond length for
specified molecule types.

Usage:

`BondLength <input.vcf> <width> <output file> <molecule names> <options>`

> `<input.vcf>`
> > input coordinate filename (must end with `.vcf`) containing indexed timesteps
> `<width>`
> > width of each bin for the distribution
> `<output file>`
> > output filename containing distribution of bond lengths
> `<molecule names>`
> > names of molecule types to calculate the distribution for

# Config utility {#Config}

This utility takes `.vcf` file containing all beads (such as `All.vcf` created by
DL_MESO `traject` utility which was modified by me)
and creates `CONFIG` file (file containing initial coordinates for a
simulation via [DL_MESO simulation
package](http://www.scd.stfc.ac.uk//research/app/ccg/software/DL_MESO/40694.aspx)).
If a `vcf` file that does not contain all beads, Config will still run, but
the generated file will not contain coordinates for all beads and thus will
not be able to be used to start a simulation rune using DL\_MESO software.

Usage:

`Config <input.vcf> <options>`

> `<input.vcf>`
> > input coordinate filename (must end with `.vcf`) containing indexed timesteps
> `<options>`
> > `-st <int>`
> > > timestep for creating the CONFIG file (if the number is higher than the
> > > number of steps, the last step is used)

# DensityAggregates {#AggDensity}

This utility calculates number bead density for aggregates of specified
size from their centre of mass. During the calculation, only the current
aggregate is taken into account, so there is no possibility of getting
'false' densities from adjacent aggregates. Therefore if some bead type is
never present in an aggregate of specified size (but is in the `.vcf` file),
its density will always be 0.

Care must be taken with beadtype names in various molecules types, because
if one beadtype appears in more molecule types, the resulting density for
that beadtype will be averaged without regard for the various types of
molecule it appears in.

Instead of true aggregate size, a number of molecules of specified name can
be used, i.e. an aggregate with 1 `A` molecule and 2 `B` molecules can be
specified with `<agg sizes>` of 3 without `-m` option or 1 if `-m A` is
used (or 2 if `-m B` is used).

Also specified molecule type(s) can be excluded via the `-x` option. This
is useful in case of several molecules sharing the same bead type.
Calculated densities take into account only name of a bead type, not in
which molecule(s) it occurs. The density from the bead type in different
molecule types will therefore be the sum of the densities from those
molecules.

Usage:

`DensityAggregates <input.vcf> <input.agg> <width> <output.rho> <agg sizes> <options>`

> `<input.vcf>`
> > input coordinate filename (must end with `.vcf`) containing indexed timesteps
> `<input.agg>`
> > input filename (must end with `.agg`) containing information about
> > aggregates
> `<width>`
> > width of each bin for the distribution
> `<output.rho>`
> > output density file (automatic ending `agg#.rho` added)
> `<agg sizes>`
> > aggregate sizes for density calculation
> `<options>`
> > `--joined`
> > > specify that the `<input.vcf>` contains aggregates with joined
> > > coordinates
> > `-n <int>`
> > > number of bins to average
> > `-st <int>`
> > > starting timestep for calculation
> > `-m <molecule type name>`
> > > instead of aggregate size, use number of molecules of specified molecule
> > > types
> > `-x <name(s)>`
> > > exclude specified molecule(s)

\todo DensityAggregates: check if only chains in one aggregate are used --
anomalies in VanDerBurgh/AddedPol/

# DensityMolecules {#DensityMolecules}

DensityMolecules works in similar way as the DensityAggregates, only instead of
aggregates, the densities are calculated for specified molecule types.

It is possible to use specified bead instead of the centre of mass for the
coordinates to calculate densities from. Care must be taken, because the order
of molecule types is taken from `vsf` rather than from `DensityMolcules`
arguments. For example: whether bead 1 will be connected with `NameA` or
`NameB` in `DensityMolecules ... NameA NameB -c 1 2` depends on molecules'
order in `vsf` file; that is if `NameA` is first in `vsf`, 1 will be
associated  with `NameA` and 2 with `NameB`, but if `NameB` is first, the
associations are reverse, regardless of the order of names in the command's
arguments. If the centre of mass should be used, `x` is given as argument. In
the above example (assuming `NameA` is first in `vsf`) if bead 1 is intended
to be used for `NameB`, but centre of mass for `NameA`, then an argument of the
form `-c x 1` must be used.

Usage:

`DensityMolecules <input.vcf> <width> <output.rho> <mol name(s)> <options>`

> `<input.vcf>`
> > input coordinate filename (must end with `.vcf`) containing indexed timesteps
> `<width>`
> > width of each bin for the distribution
> `<output.rho>`
> > output density file (automatic ending `agg#.rho` added)
> `<mol name(s)>`
> > names of molecule types to calculate density for
> `<options>`
> > `--joined`
> > > specify that the `<input.vcf>` contains aggregates with joined
> > > coordinates
> > `-n <int>`
> > > number of bins to average
> > `-c x/<int>`
> > > use specified molecule bead instead of centre of mass

# DistrAgg utility {#DistrAgg}

DistrAgg calculates number and weight average aggregation masses for each
timestep (i.e. their time evolution) as well as their average over the whole
simulation.
\latexonly
The number average and weight average aggregation mass, $\langle
M\rangle_{\mathrm{n}}$ and $\langle M\rangle_{\mathrm{w}}$ respectively, are
defined as:
\begin{equation}
  \langle M\rangle_{\mathrm{n}} = \frac{\sum_{A_{\mathrm{S}}} m_i N_i}{N} \mbox{ , \ \ and \ \ }
  \langle M\rangle_{\mathrm{w}} = \frac{\sum_{A_{\mathrm{S}}} m_i^2 N_i}{\sum_{A_{\mathrm{S}}} m_i N_i},
\end{equation}
where the sums go over all aggregate sizes. $m_i$ is the mass of an aggregate
with aggregation number $A_{\mathrm{S}}=i$, $N_i$ is the number of such
aggregates and $N$ is total number of aggregates. The equations can be
written more conveniently (for programming purposes) as sums over
aggregates themselves instead of their sizes:
\begin{equation}
  \langle M\rangle_{\mathrm{n}} = \frac{\sum_{i=1}^N m_i}{N}
  \mbox{ , \ \ and \ \ }
  \langle M\rangle_{\mathrm{w}} = \frac{\sum_{i=1}^N m_i^2}{\sum_{i=1}^N m_i}.
\end{equation} \endlatexonly

The utility also calculates number, weight and z distribution function of
aggregation numbers.
\latexonly
The distributions,
$F_{\mathrm{n}}(A_{\mathrm{S}})$, $F_{\mathrm{w}}(A_{\mathrm{S}})$, and
$F_{\mathrm{z}}(A_{\mathrm{S}})$ respecively, are defined as:
\begin{equation}
  F_{\mathrm{n}}(A_{\mathrm{S}}) =
    \frac{N_{A_\mathrm{S}}}{\sum_{A_{\mathrm{S}}} N_i} =
    \frac{N_{A_\mathrm{S}}}{N}
  \mbox{ , \ \ \ \ \ }
  F_{\mathrm{w}}(A_{\mathrm{S}}) = \frac{N_{A_\mathrm{S}}
    m_{A_\mathrm{S}}}{\sum_{A_{\mathrm{S}}} N_i m_i}
  \mbox{ , \ \ \ \ \ and \ \ \ \ \ }
  F_{\mathrm{z}}(A_{\mathrm{S}}) = \frac{N_{A_\mathrm{S}}
    m_{A_\mathrm{S}}^2}{\sum_{A_{\mathrm{S}}} N_i m_i^2},
\end{equation}
where $N_{A_{\mathrm{S}}}$ and $m_{A_{\mathrm{S}}}$ stand for the number
and mass, respectively, of aggregates with aggregation number
$A_{\mathrm{S}}$. The equations are normalized so that $\sum F_x
(A_{\mathrm{S}}) = 1$. Equations for $F_{\mathrm{w}}$ and $F_{\mathrm{z}}$
can again be transformed to contain sums over aggregates, not their sizes:
\begin{equation}
  F_{\mathrm{w}}(A_{\mathrm{S}}) = \frac{N_{A_\mathrm{S}}
    m_{A_\mathrm{S}}}{\sum_{i=1}^N m_i}
  \mbox{ \ \ \ \ \ and \ \ \ \ \ }
  F_{\mathrm{z}}(A_{\mathrm{S}}) = \frac{N_{A_\mathrm{S}}
    m_{A_\mathrm{S}}^2}{\sum_{i=1}^N m_i^2},
\end{equation} \endlatexonly

Lastly, the utility calculates distribution of volume fractions of
aggregates under the assumption that all beads have the same volume, which
is essential in dissipative particle dynamics simulation method.
\latexonly
Since in DPD the volumes of all beads are identical, volume fraction of an
aggregate with aggregation number $A_{\mathrm{s}}$ is calculated as:
\begin{equation}
\phi(A_{\mathrm{S}}) = \frac{N_{A_{\mathrm{S}}}
  n_{A_{\mathrm{S}}}}{\sum_{i=1}^N n_i},
\end{equation}
where $n_i$ is the number of beads in an aggregate with $A_{\mathrm{s}}=i$.
\endlatexonly

The definition of aggregation number\latexonly, $A_{\mathrm{S}}$,\endlatexonly
is somewhat flexible. The `-m <name>` option can be used to specify that the
aggregation number is not the number of all molecules in an aggregate, but
rather only the number molecules of the specified type(s). Two values are then
taken as the mass and aggregation number for calculations of distributions and
averages -- the 'true' value (noted as `whole agg mass` in output files) and
the '`-m` option value' (noted as `options mass`). For example, let's assume
an aggregate contains `A` and `B` molecules. Using `-m A` will count the
aggregate size as the number of `A`s, but its mass as both the mass of only `A`
molecules and as the sum of `A` and `B` molecules' masses.  The resulting
distributions will therefore be the functions of aggregation number specified
by the `-m` option, but for every distribution function, there will be two data
columns in the output file (and two sets of overall averages).

Also using the `--only <name>` uses only aggregates composed exclusively of a
specified molecule type. On the contrary, `-x <name(s)>` option discounts
aggregates containing only the specified molecules.

Lastly, the calculations can be made only for a given range of aggregation
numbers (specified by the `-m` option if present) if `-n` option is used.

The utility reads information about aggregate from input file with
[Aggregate format](\ref AggregateFile). This file can be generated using
[Aggregates utility](\ref Aggregates).

Usage:

`DistrAgg <input> <distr file> <avg file> <options>`

> `<input>`
> > input filename with information about aggregates
> `<distr file>`
> > output filename with weight and number distribution functions
> `<avg file>`
> > output filename with weight and number average aggregation number in
> > each timestep
> `<options>`
> > `-st <int>`
> > > starting timestep for calculation (does not affect calculation of
> > > time evolution)
> > `-n <int> <int>`
> > > range of aggregation numbers to calculate distributions and averages
> > > for
> > `-m <name(s)>`
> > > instead of aggregate size, use number of molecules of specified molecule
> > > types
> > `--only <molecule type name>`
> > > use only aggregates composed of specified molecule type
> > `-x <name(s)>`
> > > exclude aggregates containing only specified molecule(s)

# GyrationAggregates utility {#GyrationAggregates}

This utility calculates a gyration tensor and its eigenvalues (as roots of
the tensor's characteristic polynomial) for all aggregates. Using the
eigenvalues, various shape descriptors are determined.

It calculates radius of gyration,\latexonly $R_{\mathrm{G}}$, as:
\begin{equation}
  R_{\mathrm{G}}^2 = \lambda_x^2 + \lambda_y^2 + \lambda_z^2 \mbox{,}
\end{equation}
where $\lambda_i^2$ is the $i$-th eigenvalue of the gyration tensor (or its
$i$-th principle moment). The eigenvalues are sorted so that $\lambda_x^2
\leq \lambda_y^2 \leq \lambda_z^2$.  Then it calculates\endlatexonly the
asphericity, \latexonly $b$:
\begin{equation}
  b = \lambda_z^2 - \frac{1}{2} \left( \lambda_x^2 + \lambda_y^2 \right) =
    \frac{3}{2} \lambda_z^2 - \frac{R_{\mathrm{G}}^2}{2} \mbox{,}
\end{equation}
\endlatexonly the acylindricity, \latexonly $c$:
\begin{equation}
  c = \lambda_y^2 - \lambda_x^2
\end{equation}
\endlatexonly and the relative shape anisotropy\latexonly , $\kappa$:
\begin{equation}
  \kappa^2 = \frac{b^2 + 0.75 c^2}{R_{\mathrm{G}}^4} = \frac{3}{2}
  \frac{\lambda_x^4 + \lambda_y^4 + \lambda_z^4}{\left( \lambda_x^2 +
  \lambda_y^2 + \lambda_z^2 \right)^2}
\end{equation}
\endlatexonly

Number averages of all the properties and weight and z averages
for radius of gyration are calculated. \latexonly The number average for a quantity
$\mathcal{O}$ is defined as:
\begin{equation}
  \langle\mathcal{O}\rangle_{\mathrm{n}} = \frac{\sum_i \mathcal{O}_i}{N},
\end{equation}
where $N$ is the total number of aggregates. The weight average is:
\begin{equation}
  \langle\mathcal{O}\rangle_{\mathrm{w}} = \frac{\sum_i m_i
    \mathcal{O}_i}{\sum_i m_i},
\end{equation}
where $m_i$ is mass of an aggregate $i$.  Lastly, the z-average is:
\begin{equation}
  \langle\mathcal{O}\rangle_{\mathrm{z}} = \frac{\sum_i m^2_i
    \mathcal{O}_i}{\sum_i m^2_i}.
\end{equation}
All the sums go over individual aggregates.  Both number and weight average
according to these equations were already used in the DistrAgg utility.

\endlatexonly

Averages of the shape descriptors during the simulation (or their time
evolution) are written to an output file and overall averages are appended
to that file.

Sums dependent on aggregation number can be written to output file with
`-ps <name>` option. The aggregation number can be either the number of all
molecules in an aggregate or a number of molecules of the type specified by
`-m <name>` option.

Usage:

`GyrationAggregates <input.vcf> <input.agg> <output> <agg sizes> <options>`

> `<input.vcf>`
> > input coordinate filename (must end with `.vcf`) containing indexed timesteps
> `<input.agg>`
> > input filename (must end with `.agg`) containing information about
> > aggregates
> `<output.vcf>`
> > output filename with shape descriptors for chosen sizes throughout
> > simulation
> `<agg sizes>`
> > aggregate sizes for gyration calculation
> `<options>`
> > `--joined`
> > > specify that the `<input.vcf>` contains aggregates with joined
> > > coordinates
> > `-bt`
> > > specify bead types to be used for calculation (default is all)
> > `-m <name>`
> > > take as an aggregate size the number of `<name>` molecules in aggregates
> > > instead of the number of all molecules
> > `-ps <name>`
> > > output filename with per-size (or per-aggregation-number) averages
> > `-n <int> <int>`
> > > range of aggregation numbers to calculate distributions and averages

# GyrationMolecules utility {#GyrationMolecules}

This utility function in the same way as GyrationAggregates, but it calculates
radii of gyration for specified molecule types instead of aggregate sizes.

Right now it calculates gyration for all beads in the specified molecule types.

Usage:

`GyrationMolecules <input.vcf> <input.agg> <output> <molecule names> <options>`

> `<input.vcf>`
> > input coordinate filename (must end with `.vcf`) containing indexed timesteps
> `<input.agg>`
> > input filename (must end with `.agg`) containing information about
> > aggregates
> `<output.vcf>`
> > output filename with radii of gyration throughout simulation (automatic
> > ending #.txt)
> `<molecule names>`
> > molecule types for gyration calculation
> `<options>`
> > `-bt`
> > > specify bead types to be used for calculation (default is all)
> > `--joined`
> > > specify that the `<input.vcf>` contains joined coordinates

# JoinAggregates utility {#JoinAggregates}

This utility reads input `.vcf` and `.agg` files and removes periodic boundary
conditions from aggregates - e.i. it joins the aggregates. The distance and the
bead types for closeness check are read from the first line of `.agg` file with
contains full `Aggregates` command used to generate the file. JoinAggregates is
meant for cases, where `-j` flag was omitted in Aggregates utility.

Usage:

`Aggregates <input.vcf> <input.agg> <output.vcf> <options>`

> `<input.vcf>`
> > input coordinate filename (must end with `.vcf`) containing indexed timesteps
> `<input.agg>`
> > input filename (must end with `.agg`) containing information about
> > aggregates
> `<output.vcf>`
> > output filename (must end with `.vcf`) with joined coordinates

# JoinRuns utility {#JoinRuns}

MOST LIKELY NOT WORKING -- IT'S NOT USED.

This program is to be used if two simulation runs with different initial
seeds (that is, two simulations with different bead id numbers `.vsf`
files, but identical `FIELD` files) should be joined. Two `.vcf` files that
contain the same bead types must be provided as well as the `.vsf`
structure file for the second simulation. The output `.vcf` coordinate
files has bead ids according to the structure file of the first simulation.
The program is, however, extremely inefficient with unbonded beads, while
bonded beads are always sorted in the same way by DL_MESO simulation
software. The usefullness of such utility is confined to cases with more
then one type of unbonded beads and under those conditions the utility may
take around 1 minute per step (of the second simulation run) for system in
box of side length 40.

Usage:

`JoinRuns <1st input.vcf> <2nd input.vcf> <2nd input.vsf> <output.vcf> <type names> <options>`

>  `<1st input.vcf>`
> > input coordinate filename (must end with `.vcf`) containing indexed timesteps for the first simulation
>  `<2nd input.vcf>`
> > input coordinate filename for the second sumation in the same format as
> > the first coordinate file
>  `<2nd input.vsf>`
> > `.vsf` structure file for the second simulation (must end with `.vsf`)
>  `<output.vcf>`
> > output filename with indexed coordinates (must end with `.vcf`)
> `<type names>`
> > names of bead types to save
> `<options>`
> > `--joined`
> > > join individual molecules by removing periodic boundary conditions
> >   `-n1 <int>`
> > > number of timestep to start the first simulation from
> >   `-n2 <int>`
> > > number of timestep to start the second simulation from
> >   `-u1 <int>`
> > > leave out every `skip` steps in the first simulation
> >   `-u2 <int>`
> > > leave out every `skip` steps in the second simulation

\todo JoinRuns: base reindexing of beads in the second simulation on comparison
between the two `.vsf` files
\todo JoinRuns: implement wholy `--script` common option \todo Completely
change this - either implement `-x` option or remove function
`WriteCoorIndexed` and hard code the writing to file

# PairCorrel utility {#PairCorrel}

This utility calculates pair correlation function (pcf) between specified bead types.
All pairs of bead types (including same pair) are calculated - given `A` and `B`
types, pcf between `A-A`, `A-B` and `B-B` are calculated. The pcfs are not
normalised.

Currently, the utility cannot recognise between beads of the same type that
are in different molecule types - i.e. if bead type `A` is both in molecule
type `1` and molecule type `2`, only one pcf will be calculated regardless
of the molecule type `A` is in. Use `SelectedVcf` with `-x` option to
resolve the problem.

Usage:

`PairCorrel <input.vcf> <width> <output.pcf> <bead type(s)> <options>`

> `<input.vcf>`
> > input coordinate filename (must end with `.vcf`) containing indexed timesteps
> `<width>`
> > width of each bin for the distribution
> `<output.pcf>`
> > output file with pair correlation function(s)
> `<bead type(s)>`
> > bead type name(s) for pcf calculation
> `<options>`
> >  `-n <int>`
> > > number of bins to average
> >  `-st <int>`
> > > starting timestep for calculation

<!--
# PairCorrelPerAgg utility {#PairCorrelPerAgg}

DO NOT USE.

What would this be for? Per aggregate densities are calculated via
DensityAggregates utility.

DISREGARD THE REST OF THIS ENTRY.

PairCorrelPerAgg utility calculates pair correlation function per
aggregates - that is only beads in the same aggregate are used. If
aggregate size(s) is not specified, average pcf is calculated (that is,
regardless of aggregate size). In all probability the utility is working,
but since it is not really useful, it has never been thouroughly tested.

Usage:

`PairCorrelPerAgg <input.vcf> <input.agg> <width> <output.pcf> <bead type(s)> <options>`

> `<input.vcf>`
> > input coordinate filename (must end with `.vcf`) containing indexed timesteps
> `<input.agg>`
> > input filename (must end with `.agg`) containing information about
> `<width>`
> > width of each bin for the distribution
> `<output.pcf>`
> > output file with pair correlation function(s)
> `<bead type(s)>`
> > bead type name(s) for pcf calculation
> `<options>`
> >  `-n <int>`
> > > number of bins to average
> >  `-st <int>`
> > > starting timestep for calculation
-->

# SelectedVcf utility {#SelectedVcf}

This utility takes `.vcf` file containing either
[indexed timesteps](\ref IndexedCoorFile) and creates a new `.vcf`
coordinate file containing only beads of selected types with an option of
removing periodic boundary condition and thus joining molecules.  The otput
`.vcf` file contains indexed timesteps.

Specified molecules can be excluded which is useful when the same bead type is
shared between more molecule types.

The selected `<type names>` are printed at the beginning of the output
file.

Usage:

`SelectedVcf <input.vcf> <output.vcf> <type names> <options>`

> `<input.vcf>`
> > input coordinate filename (must end with `.vcf`) containing indexed timesteps
> `<output.vcf>`
> > output filename with indexed coordinates (must end with `.vcf`)
> `<type names>`
> > names of bead types to save
> `<options>`
> > `--join`
> > > join individual molecules by removing periodic boundary conditions
> > `-st <int>`
> > > starting timestep for calculation
> > `-sk <int>`
> > > number of steps to skip per one used
> > `-x <name(s)>`
> > > exclude specified molecule(s)

# traject utility {#traject}

This utility is from the
[DL_MESO simulation package](http://www.scd.stfc.ac.uk//research/app/ccg/software/DL_MESO/40694.aspx).
While originally it creates a `.vtf` file containing both structure and
coordinates, I have changed it to create a separate `dl_meso.vsf` structure
file and `All.vcf` coordinate file containing indexed timesteps.

There are two versions from two versions of the
[DL_MESO simulation package](http://www.scd.stfc.ac.uk//research/app/ccg/software/DL_MESO/40694.aspx),
namely versions 2.5 and 2.6.

Usage:

`traject-v2_5 <cores>` or `traject-v2_6 <cores>`

> `<cores>`
> > number of computer cores used for the simulation run (or the number of
> > `HISTORY` file)

The standard options cannot be used with this utility.> }}}

# TransformVsf utility {#TransformVsf}

This utility takes `.vsf` structure file and DL_MESO input file `FIELD` and
transforms them into a different `.vsf` structure file that is better
suited for visualisation using VMD software.

Usage:

`TransformVsf <output.vsf> <options>`

> `<output.vsf>`
> > output structure file that must end with `.vsf`

## Format of output structure file {#TransformVsf-output}

Every atom line in the generated structure file contains bead's index number,
mass, charge and name. Atom lines for beads in molecules also contain
molecule's id number and the name of the type of molecule. The bond section of
`output.vsf` lists all bonds one by one (i.e. no chains of bonds in the format
`<id1>:: <id2>` are used). Information about which bonds belong to which
molecule is provided as a comment. The file has the following format:
> `atom default name <name> mass <m> charge <q>`
>
> `...`
>
> `atom <id> name <name> mass <m> charge <q>`
>
> `...`
>
> `atom <id> name <name> mass <m> charge <q> segid <name> resid <id>`
>
> `...`
> 
> `# resid <id>`
>
> `<bonded bead id1>: <bonded bead id2>`
>
> `...`

For VMD atom selection:
> `segid <name>`
> > selects all molecules with given name(s)
> `resid <id>`
> > selects molecule(s) with given index number(s)
> `charge <q>`
> > selects all beads with given charge(s) (double quotes are required for
> > negative charge)
> `mass <m>`
> > selects all beads with given mass(es)
