#!/bin/bash -e

Dir=`pwd`

for pH in 08.00 09.00 10.00
do
    mkdir -p occ_pH${pH}
    cp ../../03_CpHMetaD/pH${pH}/U1mer_???.occ occ_pH${pH}
done
