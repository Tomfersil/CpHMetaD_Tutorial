#!/bin/bash -e

Dir=`pwd`
sys=U1mer
begin_time=2

for pH in 08.00 09.00 10.00
do

	gmx trjcat -f ${sys}_pH${pH}.xtc -o tmp_pH${pH}.xtc -b ${begin_time}

done

gmx trjcat -f tmp_pH* -o ${sys}_tot.xtc -cat

