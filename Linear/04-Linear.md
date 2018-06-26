# Utilities for linear chains {#LinearChains}

This section provides information about utilities with calculations that
are sensible to do only on linear polymer chains. No check whether the
molecules are linear is done.

# EndToEnd utility {#EndToEnd}

This utility calculates end-to-end distance of specified molecules.
End-to-end distance makes sense only for linear chains, therefore it is
assumed that the provided molecule names are linear chains. No check is
performed. The distance is calculated between the first and the last bead
of the molecule; that is, between the first and the last bead in the
`FIELD` entry for the given molecule. Also the use of joined coordinates
(that is, without periodic boundary condition) is required, because the
utility does not remove periodic boundary conditions.

The output is a file containing average end-to-end distance for every
molecule type for each timestep.

Usage:

`EndToEnd <input.vcf>  <output.vcf> <molecule names> <options>`

> `<input.vcf>`
> > input coordinate filename (must end with `.vcf`) containing either
> > ordered or indexed timesteps (with joined coordinates)
> `<output.vcf>`
> > output filename with indexed coordinates (must end with `.vcf`)
> `<molecule names>`
> > names of molecule types (linear chains) to use
