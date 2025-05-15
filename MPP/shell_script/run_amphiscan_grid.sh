#!/bin/bash
#

start=`date +%s`

if [ -d results ]; then
	rm -r results
fi

rm rmsd_values.csv
mkdir results

ls input_pdbs/* > files
sed -i 's=input_pdbs/==g' files
sed -i 's=.pdb==g' files
mapfile -t array <files

for element in ${array[@]}; do
	# Run AmphiScan0_2.py with the current z_increment and pert_num
	python311 AmphiScan0_2.py -input_pdb ${element}.pdb

	# Calculate RMSD with the current z_increment and pert_num
	python311 rmsd.py -protein ${element}.pdb
done
python311 rmsd_plot.py -pert_num 15 -z_increment 1
end=`date +%s`
runtime=$((end-start))
hours=$((runtime / 3600)); minutes=$(( (runtime % 3600) / 60 )); seconds=$(( (runtime % 3600) % 60 ));
echo "Runtime: $hours:$minutes:$seconds (hh:mm:ss)"