#!/bin/bash

# 1) surround part of the diblock with a layer of monomer beads
AddToSystem in.vcf 1.v{s,c}f -f FIELD-1 -ld 1.5 -hd 1.501 -bt B -i in.vsf -cy 5.5 10 -cz 2 6.2
# 2) switch part of the surrounding beads for molecules placed randomly in a bigger box
AddToSystem 1.vcf 2.v{s,c}f -f FIELD-2 -i 1.vsf --switch -b 10 15 10
# 3) add monomer beads to the lower part of the box
AddToSystem 2.vcf 3.v{s,c}f -f FIELD-3 -i 2.vsf -cy 0 5
# 4) add the system from 1) to 3) while enlarging the box and placing 3) in the middle of the new box
AddToSystem 3.vcf 4.v{s,c}f -vtf 1.v{s,c}f -b 20 15 10 --centre -i 3.vsf
