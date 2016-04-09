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

\todo Implement possibility to choose timestep number for creating `CONFIG`
file.

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

\todo Somehow avoid the need to use the special optional bond file, where
the ids of beads must strictly adhere to FIELD. Possibly require use of a
vcf file in conjunction with bond file

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

\todo Add the possibility to save only certain bead types to output vcf
file with joined coordinates.

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

\todo Look into the number averages.

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

\todo Completely rewrite - especially remove requirement for number of
lines on the first line of input file
