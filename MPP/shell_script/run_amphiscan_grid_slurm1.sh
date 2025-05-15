#!/bin/bash

start=`date +%s`
cd /home/cadeniran/scratch/amphiscan-pep
current_dir=$PWD

# Define grid search parameters for z_increment and pert_num
z_increments=(1.0)
#z_increments=(0.9 0.8 0.7)
pert_nums=(10 15 20 25 30 35 40 45 50)

# Iterate over each PDB file and submit a SLURM job for each one
for z_increment in "${z_increments[@]}";
do
	for pert_num in "${pert_nums[@]}";
	do
		if [ ! -d results_grid/${z_increment}_${pert_num} ];
		then
			mkdir -p results_grid/${z_increment}_${pert_num}
			cp -R Amphiscan_rot.py modules.py input_pdbs rmsd_plot.py testrmsd.py testembedres.py testembedres_AT.py membrane_thickness -t results_grid/${z_increment}_${pert_num}
			sed -e 's|\$zinc|'"$z_increment"'|g' -e 's|\$pert|'"$pert_num"'|g' t > results_grid/${z_increment}_${pert_num}/run_amphiscan0.2.sh
			cd results_grid/${z_increment}_${pert_num}
			nohup bash run_amphiscan0.2.sh &
			echo "Submitted: ${z_increment}_${pert_num}"
			sleep 2
			cd "$current_dir"
		fi
	done
done
