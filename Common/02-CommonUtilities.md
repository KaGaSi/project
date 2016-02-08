# Common utilities {#Common}

Utilies that are not specific to any given system and are used for all
simulations.

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
> > `-i <name>`
> > > use custom `.vsf` structure file instead of the default `dl_meso.vsf`
> > > (must end with `.vsf`)
> > `-b <name>`
> > > file containing bond alternatives to `FIELD`
> > `-j`
> > > join individual molecules by removing periodic boundary conditions
> > `-v`
> > > verbose output providing information about the system
> > `-V`
> > > more detailed verbose output (also prints comments from `.vcf` file
> > > at the start of every timestep)
> > `-h`
> > > print help and exit

\todo Implement `<skip>` and `<start>`

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
> `<options>`
> > `-i <name>`
> > > use custom `.vsf` structure file instead of the default `dl_meso.vsf`
> > > (must end with `.vsf`)
> > `-b <name>`
> > > file containing bond alternatives to `FIELD`
> > `-v`
> > > verbose output providing information about the system
> > `-V`
> > > more detailed verbose output (also prints comments from `.vcf` file
> > > at the start of every timestep)
> > `-h`
> > > print help and exit

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
> `<options>`
> > `-i <name>`
> > > use custom `.vsf` structure file instead of the default `dl_meso.vsf`
> > > (must end with `.vsf`)
> > `-b <name>`
> > > file containing bond alternatives to `FIELD`
> > `-v`
> > > verbose output providing information about the system
> > `-h`
> > > print help and exit

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
> `<output.agg>`
> > output filename (must end with `.agg`) containing information about
> > aggregates
> `<distance>`
> > minimum distance for two beads to be in contact (constituting one
> > contact pair)
> `<contacts>`
> > minimum number of contact pairs to consider two molecules to be in one
> > aggregate
> `<type names>`
> > names of bead types to use for calculating contact pairs
> `<options>`
> > `-i <name>`
> > > use custom `.vsf` structure file instead of the default `dl_meso.vsf`
> > > (must end with `.vsf`)
> > `-b <name>`
> > > file containing bond alternatives to `FIELD`
> > `-j <joined.vcf>`
> > > filename for coordinates of joined aggregates (must end with `.vcf`)
> > `-v`
> > > verbose output providing information about the system
> > `-V`
> > > more detailed verbose output (also prints comments from `.vcf` file
> > > at the start of every timestep)
> > `-h`
> > > print help and exit

\todo Add the possibility to save only certain bead types to output vcf
file with joined coordinates.
