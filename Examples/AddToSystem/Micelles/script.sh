#!/bin/bash

bin=~/AnalysisTools/build/bin

# 1) create the first molecule
#      -cx, -cy, and -cz constraints ensure the first bead is in the middle of the box with size 5x5x5
$bin/AddToSystem -- 1.v{s,c}f -b 5 5 5 -cx 2.49 2.51 -cy 2.49 2.51 -cz 2.49 2.51 -f FIELD-1
# 2) add the remaining 19 molecules and generate FIELD for the micelle
#      -ld, -hd, and -bt ensure evry molecule's first bead is close to that of the molecule from 1)
#      -f_out creates a FIELD file containing the micelle as a single molecule
$bin/AddToSystem 1.vcf 2.v{s,c}f -i 1.vsf -ld 0.5 -hd 0.51 -bt A -f FIELD-2 -f_out micelle.FIELD 1
# 3) populate a box with 5 micelles
#     first, place randomly the first micelle into a box of size 20 20 5
#     second, place into the box one by one another 4 micelles, ensuring the
#       micelles do not overlap by using -ld and -bt options
#     -cx and -cy constraints ensure the micelle does not stick out of the box
out=b
$bin/AddToSystem -- $out.v{s,c}f -b 20 20 5 -gc -f micelle.FIELD -cx 2.5 17.5 -cy 2.5 17.5
for (( i=0 ; i<4 ; i++ )); do
  in=a
  mv $out.vsf $in.vsf
  mv $out.vcf $in.vcf
  $bin/AddToSystem $in.vcf $out.v{s,c}f -i $in.vsf -ld 4 -bt A B C -gc -f micelle.FIELD -cx 2.5 17.5 -cy 2.5 17.5
  rm $in.v{s,c}f
  sleep 1
done
mv $out.vsf 3.vsf
mv $out.vcf 3.vcf
