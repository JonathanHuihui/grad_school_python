#!/bin/bash
# Take every (i+10)th frame from the trajectory and compute the N by N matrix for that frame
num_frames=100000

topfile=$( ls -v *.top )
org=${topfile%.*}
startRow=$(( $( awk '/FLAG ATOM_NAME/{ print NR }' ${topfile} ) + 2 ))
stopRow=$(( $( awk '/FLAG CHARGE/{ print NR }' ${topfile} ) - 1 ))
NRESI=$( sed -n "$startRow","$stopRow"p ${topfile} | grep -o "CA" | wc -l )


for (( iter=0; iter<=100000; iter+=10 )); 
do
	echo trajin /home/amber/resurrected_trx/aeca/aeca_1000ns.mdcrd ${iter} $((${iter} + 1)) 1 >> temp.in
	echo "" >> temp.in
	echo matrix out aeca_distance_${iter}.dat name aeca_distance_${iter} byres :1-${NRESI} dist >> temp.in
	cpptraj -i temp.in -p ${topfile} 
	rm temp.in
done

tar -cvf ${org}.tar *.dat
rm *.dat

