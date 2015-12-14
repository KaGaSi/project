Program TransformVsf {#TransformVsf}
=====

This program takes `.vsf` structure file and DL_MESO input file `FIELD` and
transforms them into a different `.vsf` structure file that is well suited
for visualisation using VMD software. The output file contains name, mass
and charge of every bead. In case of a bead in a molecule, it also contains
its name and id.

Required format of input .vsf structure file {#TransformVsf-Require}
=====

The program is designed with file `dl_meso.vsf` in mind, which is generated
by the `traject` program provided in DL_MESO software (and modified by me).
`TransformVsf` was tested only against files generated by `traject`, but
other `.vsf` files should work fine, if formatted according to the
following guidelines.

The first line specifies default bead type which means all atom lines for
beads of this type are unnecessary (provided those beads are not in a
molecule). All atom lines specify VDW radius and atom name. If an atom is
in a residue, its residue number is appended to the atom line.
> `atom default radius 1.000000 name <name>`
>
> `atom <id> radius 1.000000 name <name> resid <id>`
>
> `...`

Only the bead number and name are read, so both VDW radius and residue
number are not strictly necessary.  Short version of `atom` and `name`
keywords (`a` and `n` respectively) can be used.  Other keywords can be
included, since `TransformVsf` will disregard them. No comments are
allowed.

Bond lines of `.vsf` are not read and are therefore irrelevant to
`TransformVsf`.

Format of output .vsf structure file {#TransformVsf-output}
=====

Every atom line in the generated structure file contains bead's index
number, mass, charge and name. Atom lines for beads in molecules also
contain molecule's id number and the name of the type of molecule. The
bond section of `output.vsf` lists all bonds one by one (i.e. no chains of
bonds in the format `<id1>:: <id2>` are used). The file has the
following format:
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
> `<bonded bead id1>: <bonded bead id2>`
>
> `...`

For VMD atom selection:
* `segid <name>` selects all molecules with given name(s)
* `resid <id>` selects molecule(s) with given index number(s)
* `charge <q>` selects all beads with given charge(s) (double quotes are
  required for negative charge)
* `mass <m>` selects all beads with given mass(es)