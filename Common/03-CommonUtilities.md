# Common utilities {#Common}

Utilies that are not specific to any given system and are used for all
simulations.

# traject utility {#traject}

This utility is from the
[DL_MESO simulation package](http://www.scd.stfc.ac.uk//research/app/ccg/software/DL_MESO/40694.aspx).
While originally it creates a `.vtf` file containing both structure and
coordinates, I have changed it to create a separate `dl_meso.vsf` structure
file and `All.vcf` coordinate file containing ordered timesteps.

Usage:

`traject <cores>`

> `<cores>`
> > number of computer cores used for the simulation run (or the number of
> > `HISTORY` file)

The standard options cannot be used with this utility.

# SelectedVcf utility {#SelectedVcf}

This utility takes `.vcf` file containing either
[ordered timesteps](\ref OrderedCoorFile) (such as `All.vcf` created by
DL_MESO `traject` utility which was modified by me) or
[indexed timesteps](\ref IndexedCoorFile)
and creates a new `.vcf` coordinate file containing only beads
of selected types with an option of removing periodic boundary condition
and thus joining molecules.  The otput `.vcf` file therefore contains
indexed timesteps.

Usage:

`SelectedVcf <input.vcf> <start> <skip> <output.vcf> <type names> <options>`

> `<input.vcf>`
> > input coordinate filename (must end with `.vcf`) containing either
> > ordered or indexed timesteps
> `<start>`
> > number of timestep to start from
> `<skip>`
> > leave out every `skip` steps
> `<output.vcf>`
> > output filename with indexed coordinates (must end with `.vcf`)
> `<type names>`
> > names of bead types to save
> `<options>`
> > `-j`
> > > join individual molecules by removing periodic boundary conditions

# Config utility {#Config}

This utility takes `.vcf` file containing either
[ordered timesteps](\ref OrderedCoorFile) (such as `All.vcf` created by
DL_MESO `traject` utility which was modified by me) or
[indexed timesteps](\ref IndexedCoorFile)
and creates `CONFIG` file (file containing initial coordinates for a
simulation via [DL_MESO simulation
package](http://www.scd.stfc.ac.uk//research/app/ccg/software/DL_MESO/40694.aspx)).

Usage:

`Config <input.vcf> <options>`

> `<input.vcf>`
> > input coordinate filename (must end with `.vcf`) containing either
> > ordered or indexed timesteps

# TransformVsf utility {#TransformVsf}

This utility takes `.vsf` structure file and DL_MESO input file `FIELD` and
transforms them into a different `.vsf` structure file that is well suited
for visualisation using VMD software.

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
molecule is provided as has comment. The file has the following format:
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

\todo Optional bond file format: somehow avoid the need to use the special
optional bond file, where the ids of beads must strictly adhere to FIELD.
Possibly require use of a vcf file in conjunction with bond file

# BondLength utility {#BondLength}

BondLength utility calculates normalized distribution of bond length for
specified molecule types.

Usage:

`BondLength <input.vcf> <output file> <width> <molecule names> <options>`

> `<input.vcf>`
> > input coordinate filename (must end with `.vcf`) containing either
> > ordered or indexed timesteps
> `<output file>`
> > output filename containing distribution of bond lengths
> `<width>`
> > width of each bin for the distribution
> `<molecule names>`
> > names of molecule types to calculate the distribution for

# Aggregates utility {#Aggregates}

This utility determines which molecules belong to which aggregates
according to a simple criterion: two molecules belong to the same aggregate
if they at least a specified number of contact pairs. A contact pair is a
pair of two beads belonging to different molecules which are closer than
certain distance. Both the distance and the number of needed contact pairs
are arguments of the command.

Usage:

`Aggregates <input.vcf> <distance> <contacts> <output.agg> <type names> <options>`

> `<input.vcf>`
> > input coordinate filename (must end with `.vcf`) containing either
> > ordered or indexed timesteps
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
> > `-j <joined.vcf>`
> > > filename for coordinates of joined aggregates (must end with `.vcf`)

The NotSameBeads variant of the Aggregate utility works in exactly the same,
but does not calculates contacts between beads of the same type, i.e. if bead
types `A` and `B` are provided, Aggregates will calculate contact pairs `A-B`,
`A-A` and `B-B` (provided the beads are in different molecules), while
Aggregates-NotSameBeads will calculate only `A-B` contact pair. Therefore at
least two bead types must be provided for `<type names>` argument.

# JoinAggregates utility {#JoinAggregates}

This utility reads input `.vcf` and `.agg` files and removes periodic
boundary conditions from aggregates - e.i. it joins the aggregates. The
distance and the bead types for closeness check are read from the first
line of `.agg` file with contains full Aggregates command used to generate
the file. JoinAggregates is meant for cases, where `-j` flag was omitted
in Aggregates utility.

Usage:

`Aggregates <input.vcf> <input.agg> <output.vcf> <options>`

> `<input.vcf>`
> > input coordinate filename (must end with `.vcf`) containing either
> > ordered or indexed timesteps
> `<input.agg>`
> > input filename (must end with `.agg`) containing information about
> > aggregates
> `<output.vcf>`
> > output filename (must end with `.vcf`) with joined coordinates

# DistrAgg utility {#DistrAgg}

DistrAgg calculates number and weight average aggregation numbers
for each timestep.
\latexonly
The number average aggregation number, $\langle A_{\mathrm{s}} \rangle_n$
is defined as:

\begin{equation}
\langle A_{\mathrm{s}} \rangle_n = \frac{\sum_i m_i}{N} \mbox{,}
\end{equation}

where $m_i$ is weight (aggregation number) of aggregate $i$ and $N$ is
total number of aggregates. The weight average aggregation number, $\langle
A_{\mathrm{s}} \rangle_w$ is then defined as:

\begin{equation}
\langle A_{\mathrm{s}} \rangle_w = \frac{\sum_i m_i^2}{\sum_i m_i} \mbox{.}
\end{equation}
\endlatexonly

It also calculates overall number and weight distribution function.
\latexonly
The number distribution function, $F_n (A_{\mathrm{s}})$ is defined as:

\begin{equation}
F_n (A_{\mathrm{s}}) = \frac{N_{A_{\mathrm{s}}}}{\sum_i N_i} \mbox{,}
\end{equation}

where $N_i$ is the number of aggregates with aggregation number
$A_{\mathrm{s}} = i$.  The weight distribution function, $F_w
(A_{\mathrm{s}})$ is then defined as:

\begin{equation}
F_w(A_{\mathrm{s}}) = \frac{m_{A_{\mathrm{s}}} N_{A_{\mathrm{s}}}}{\sum_i m_i N_i} \mbox{,}
\end{equation}

where $m_{A_{\mathrm{s}} }$ is the weight, that is the number of aggregates
with the given aggregation number.
\endlatexonly

Lastly, the utility calculates volume fractions of all aggregates, where it
(for now) assumes that all beads have reduced mass of 1.
\latexonly
Volume fraction of an aggregate with aggregation number $A_{\mathrm{s}}$ is
defined as:

\begin{equation}
\phi(A_{\mathrm{s}}) = \frac{m_{A_{\mathrm{s}}} N_{A_{\mathrm{s}}}}{\sum_i m_i N_i} \mbox{,}
\end{equation}

where $m_i$ is the actual mass of an aggregate -- it equals to aggregate's
volume assuming all beads have a unit mass.
\endlatexonly

The utility reads information about aggregate from input file with
[Aggregate format](\ref AggregateFile)

This file can be generated using [Aggregates utility](\ref Aggregates).

Usage:

`DistrAgg <input> <output distr file> <output avg file> <options>`

> `<input>`
> > input filename with information about aggregates
> `<output distr file>`
> > output filename with weight and number distribution functions
> `<output avg file>`
> > output filename with weight and number average aggregation number in
> > each timestep
> `<options>`
> > `-n <int>`
> > > starting timestep for calculation

\todo DistrAgg: look into volume fractions with beads of arbitrary (and
different) masses.

# DensityAggregates {#AggDensity}

This utility calculates number bead density for aggregates of specified
size from their center of mass. During the calculation, only the current
aggregate is taken into account, so there is no possibility of getting
'false' densities from adjacent aggregates. Therefore if some bead type is
never present in an aggregate of specified size, its density will always be
0.

Usage:

`DensityAggregates <input.vcf> <input.agg> <width> <output.rho> <agg sizes> <options>`

> `<input.vcf>`
> > input coordinate filename (must end with `.vcf`) containing either
> > ordered or indexed timesteps
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
> > `-j`
> > > specify that the `<input.vcf>` contains aggregates with joined
> > > coordinates

\todo DensityAggregates: check if only chains in one aggregate are used --
anomalies in VanDerBurgh/AddedPol/

# DensityMolecules {#DensityMolecules}

DensityMolecules works in similar way as the DensityAggregates, only instead of aggregates,
the densities are calculated for specified molecule types.  Care must be taken
with beadtype names in various molecules types, because if one beadtype appears
in more molecule types, the resulting density for that beadtype will be
averaged without regard for the various types of molecule it appears in.

Usage:

`DensityMolecules <input.vcf> <input.agg> <width> <output.rho> <agg sizes> <options>`

> `<input.vcf>`
> > input coordinate filename (must end with `.vcf`) containing either
> > ordered or indexed timesteps
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
> > `-j`
> > > specify that the `<input.vcf>` contains aggregates with joined
> > > coordinates

# GyrationAggregates utility {#GyrationAggregates}

This utility calculates a gyration tensor and its eigenvalues (using Jacobi
transformations) for aggregates of given size. It then determines their
radius of gyration. It saves average radii of gyration during the
simulation to an output file and prints overall average.

Right now it calculates gyration for all beads in the aggregates.

Usage:

`GyrationAggregates <input.vcf> <input.agg> <output> <agg sizes> <options>`

> `<input.vcf>`
> > input coordinate filename (must end with `.vcf`) containing either
> > ordered or indexed timesteps
> `<input.agg>`
> > input filename (must end with `.agg`) containing information about
> > aggregates
> `<output.vcf>`
> > output filename with radii of gyration throughout simulation (automatic
> > ending #.txt)
> `<agg sizes>`
> > aggregate sizes for gyration calculation
> `<options>`
> > `-j`
> > > specify that the `<input.vcf>` contains aggregates with joined
> > > coordinates

\todo GyrationAggregates: implement using only selected bead types for
calculation

\todo GyrationAggregates: implement option for using only specific molecule
names for finding out the aggregation number (i.e. make
OneTime/GyragionAggregates-NotHomopol more general)

\todo GyrationAggregates: understand `jacobi` function

# GyrationMolecules utility {#GyrationMolecules}

This utility function in the same way as GyrationAggregates, but it
calculates radii of gyration for specified molecule names instead of
aggregate sizes.

Right now it calculates gyration for all beads in the aggregates.

Usage:

`GyrationMolecules <input.vcf> <input.agg> <output> <molecule names> <options>`

> `<input.vcf>`
> > input coordinate filename (must end with `.vcf`) containing either
> > ordered or indexed timesteps
> `<input.agg>`
> > input filename (must end with `.agg`) containing information about
> > aggregates
> `<output.vcf>`
> > output filename with radii of gyration throughout simulation (automatic
> > ending #.txt)
> `<molecule names>`
> > molecule types for gyration calculation
> `<options>`
> > `-j`
> > > specify that the `<input.vcf>` contains joined coordinates

\todo GyrationMolecules: implement using only selected bead types for
calculation

\todo GyrationMolecules: understand `jacobi` function

# JoinRuns utility {#JoinRuns}

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
> > input coordinate filename (must end with `.vcf`) containing either
> > ordered or indexed timesteps for the first simulation
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
> > `-j`
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

# Average utility {#Average}

Utility calculating average values with standard deviation and
autocorrelation time from values contained in a text file. The first line
of the file has to contain the number of data lines and no comments are
allowed.

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

\todo Average: completely rewrite - especially remove requirement for number of
lines on the first line of input file
