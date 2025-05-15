#!/bin/bash

dir=$PWD

if [ -d results_grid ]; then
	rm -r results_grid
fi

mkdir -vp results_grid

# Define array
myArray=("30 25 10"
"30 25 5"
"30 20 10"
"30 20 5"
"30 15 5"
"30 10 5"
"25 15 5"
"25 10 5"
"20 15 10"
"20 10 5"
"20 5  1"
"10 5  1"
"5 3 1"
"3 2 1")

# Iterate over each PDB file and submit a SLURM job for each one
while read -r a b c
do
	if [ ! -d results_grid/${a}_${b}_${c} ]; then
		mkdir -vp results_grid/${a}_${b}_${c}
		cp -R {modules.py,rmsd.py,rmsd_plot.py,input_pdbs,run_amphiscan_grid.sh} results_grid/${a}_${b}_${c}
		sed -e 's|\$rot1|'"$a"'|g' -e 's|\$rot2|'"$b"'|g' -e 's|\$rot3|'"$c"'|g' t > results_grid/${a}_${b}_${c}/AmphiScan0_2.py
		cd results_grid/${a}_${b}_${c}
	fi
	nohup ./run_amphiscan_grid.sh &
	sleep 3
	cd $dir

done < <(printf '%s\n' "${myArray[@]}")

