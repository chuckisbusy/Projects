#!/bin/bash

start=`date +%s`
cd /home/cadeniran/scratch/amphiscan-pep
current_dir=$PWD

# Define grid search parameters for z_increment and pert_num
z_increments=(0.1)
#z_increments=(0.9 0.8 0.7)
pert_nums=(170 175 180 185 190 195 200 205 210 215 220 225 230 235 240 245 250 255 260 265 270 275 280 285 290 295 300 305 310 315 320 325)

# Iterate over each PDB file and submit a SLURM job for each one
for z_increment in "${z_increments[@]}";
do
	for pert_num in "${pert_nums[@]}";
	do
		if [ ! -d results_grid/${z_increment}_${pert_num} ];
		then
			mkdir -p results_grid/${z_increment}_${pert_num}
			cp -R Amphiscan_perturb.py modules.py input_pdbs rmsd_plot.py testrmsd.py testembedres.py testembedres_AT.py membrane_thickness -t results_grid/${z_increment}_${pert_num}
			sed -e 's|\$zinc|'"$z_increment"'|g' -e 's|\$pert|'"$pert_num"'|g' t > results_grid/${z_increment}_${pert_num}/run_amphiscan0.2.sh
			cd results_grid/${z_increment}_${pert_num}
			nohup bash run_amphiscan0.2.sh &
			echo "Submitted: ${z_increment}_${pert_num}"
			sleep 2
			cd "$current_dir"
		fi
	done
done
