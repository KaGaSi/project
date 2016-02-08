# Utilities for nanoparticles {#Nanoparticles}

This section provides information about utilities used to handle
simulations with nanoparticles.

# C60 utility {#C60}

This program creates a nanoparticle of a given radius. The inner part is
formed by a C60 fullerene with diameter of 0.5 of the nanoparticle
radius. The shell is formed by a C60 fullerene with added beads in the
middle of every C60 face. There 153 beads in the nanoparticle.

Bonds are created between all beads, that is 153 * 152 / 2 = 11628 bonds
in all.

Two files are created: `C60.vtf` to check how the nanoparticle looks and
`Nano-FIELD` with coordinates and bonds in the format required by the
`FIELD` file in
[DL_MESO simulation package](http://www.scd.stfc.ac.uk//research/app/ccg/software/DL_MESO/40694.aspx)

Usage:

`C60 <radius>`

> `<radius>`
> > required radius of nanoparticle (in reduced units)
