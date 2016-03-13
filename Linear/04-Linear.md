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

# PersistenceLength utility {#PersistenceLength}

This utility calculates persistence length of specified molecules.
It is
assumed that the provided molecules are linear chains, but no check is
performed.
Also the use of joined coordinates
(that is, without periodic boundary condition) is required, because the
utility does not remove periodic boundary conditions.

The calculation of the persistence length is based on the projection of
angles between bonds vectors (see e.g.
[this paper](http://pubs.acs.org/doi/full/10.1021/ma012052u)).
\latexonly
The following formula for the persistence length, $l_{\mathrm{t}}$ is used:

\begin{equation}
  l_{\mathrm{P}} = \langle b \rangle \sum_{i=0}^{i=N_b} \langle \cos
  \theta_i \rangle \mbox{,}
\end{equation}

where $\langle b \rangle$ is the average bond length in a molecule,
$\langle \theta_i \rangle$ is the average angle between two bond vectors
separated by $i$ bonds. $N_b$ is the number of bonds in the given molecule.
\endlatexonly

The output is a file containing average persistence length for every
molecule type for each timestep.

Usage:

`PersistenceLength <input.vcf> <output.vcf> <molecule names> <options>`

> `<input.vcf>`
> > input coordinate filename (must end with `.vcf`) containing either
> > ordered or indexed timesteps (with joined coordinates)
> `<output.vcf>`
> > output filename with indexed coordinates (must end with `.vcf`)
> `<molecule names>`
> > names of molecule types (linear chains) to use
