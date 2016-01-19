Aggragates utility {#Aggregates}
=====

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
