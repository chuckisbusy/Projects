#!/bin/bash
#

start=`date +%s`
rm -r results

ls input_pdbs/* > files
sed -i 's=input_pdbs/==g' files
sed -i 's=.pdb==g' files

mapfile -t array <files
for element in "${array[@]}"; do
	python311 AmphiScan0_2.py -input_pdb ${element}.pdb -pert_num 10 -z_increment 0.25
	python311 rmsd.py -protein ${element}.pdb
	python311 rmsd_plot.py -pert_num 10 -z_increment 0.25
done

end=`date +%s`
runtime=$((end-start))
hours=$((runtime / 3600)); minutes=$(( (runtime % 3600) / 60 )); seconds=$(( (runtime % 3600) % 60 ));
echo "Runtime: $hours:$minutes:$seconds (hh:mm:ss)"
