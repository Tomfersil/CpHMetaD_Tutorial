Dir=`pwd`
sys=U1mer
tpr="../../03_CpHMetaD/pH00.00/${sys}_001.tpr"
ndx=index.ndx
equil=000000



for pH in -01.00 00.00 01.00  
do

echo $pH
rm -f ${sys}_pH${pH}.xtc

gmx trjcat -f ../../03_CpHMetaD/pH${pH}/${sys}_*.xtc -o tmp_${pH}.xtc

gmx trjconv -f tmp_${pH}.xtc -n $ndx -s $tpr -o whole_${pH}.xtc -pbc whole -b $equil <<EOF
RNA
EOF

gmx trjconv -f whole_${pH}.xtc -n $ndx -s $tpr -o ${sys}_pH${pH}.xtc -center -pbc mol <<EOF
RNA
RNA
EOF

rm -f whole_${pH}.xtc tmp_${pH}*.xtc

cat ../../03_CpHMetaD/pH${pH}/HILLS_{001..010} > HILLS_pH${pH}

rm -f whole_${pH}_part?.xtc tmp_${pH}.xtc
    
done
