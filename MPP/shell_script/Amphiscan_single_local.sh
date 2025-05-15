#!/bin/bash
start=`date +%s`

if [ -d results ]; then
rm -r results
fi

mkdir results

mapfile -t array <files

for element in ${array[@]}; do
# Run AmphiScan0_2.py with the current z_increment and pert_num
/data/sbgrid/x86_64-linux/pyrosetta/2022.12/bin/python AmphiScan_perturb.py -input_pdb ${element}.pdb -z_increment 0.5 -pert_num 3
# Calculate RMSD with the current z_increment and pert_num
/data/sbgrid/x86_64-linux/pyrosetta/2022.12/bin/python testrmsd.py -protein ${element}
done

end=`date +%s`
runtime=$((end-start))
hours=$((runtime / 3600)); minutes=$(( (runtime % 3600) / 60  )); seconds=$(( (runtime % 3600 ) % 60  ));
echo "$SLURM_JOB_NAME Runtime: $hours:$minutes:$seconds (hh:mm:ss)" >> jobtime.txt
