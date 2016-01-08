SelectedVcf utility {#SelectedVcf}
=====

This utility takes `.vcf` file containing either
[ordered timesteps](\ref OrderedCoorFile) (such as `All.vcf` created by
DL_MESO `traject` utility which was modified by me) or
[indexed timesteps](\ref IndexedCoorFile)
and creates a new `.vcf` coordinate file containing only beads
of selected types.  The otput `.vcf` file therefore contains indexed timesteps.

Usage:

`SelectedVcf <input.vcf> <start> <skip> <output.vcf> <type names> <options>`

> `<input.vcf>`
> > input coordinate filename (must end with `.vcf` and contain ordered
> > timesteps)
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
> > `-v`
> > > verbose output providing information about the system
> > `-V`
> > > more detailed verbose output (also prints comments from `.vcf` file
> > > at the start of every timestep)
> > `-h`
> > > print help and exit

\todo Implement `<skip>` and `<start>`
