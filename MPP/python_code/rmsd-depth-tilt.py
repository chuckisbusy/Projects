from pyrosetta import * 
from modules import analyze_y
from modules import sns_heatmap
from modules import AH_CA_rmsd
from modules import HelixTools
import argparse
import re
import numpy as np
import pandas as pd
import os

init()

subdir=[]
dir='/data/storage/cadeniran/mpp/memscan-local/results_grid/'
subdir = [ item for item in os.listdir(dir) if os.path.isdir(os.path.join(dir, item)) ]
subdir.sort()

proteinlist=[]
path='/data/storage/cadeniran/mpp/memscan-local/input_pdbs/'
proteinlist = [f for f in os.listdir(path) if os.path.isfile(os.path.join(path, f))]
proteinlist.sort()

depthlist = []
rmsdlist=[]
tiltlist=[]

for i in subdir:
    os.chdir(os.path.join(dir,i))
    for j in proteinlist:
        protein = j
        protein_tag = protein.split(sep='_')[0]
        multiple_tag = protein.split(sep='_')[1]

        native = pose_from_pdb('input_pdbs/{}_{}_renum.pdb'.format(protein_tag,multiple_tag))
        best = pose_from_pdb('results/{}_{}/output_pdbs/{}_{}_renum_best_pose_overall.pdb'.format(protein_tag,multiple_tag,protein_tag,multiple_tag))

        ### DEPTH CALC
        regexp = re.compile(r'The best depth.*?([0-9.-]+)')
        with open('results/{}_{}/txt/{}_{}_renum_best_scores.txt'.format(protein_tag,multiple_tag,protein_tag,multiple_tag)) as f:
            for line in f:
                match = regexp.match(line)
                if match:
                    #print(match.group(1))
                    depthlist.append(match.group(1))


        ### RMSD CALC
        rmsd = AH_CA_rmsd(native, best)
        rmsdlist.append(rmsd)


        ### TILT CALC
        v1 = HelixTools().calculate_screw_axis(best)
        bestdegree = HelixTools().calc_angle(v1, 'z')

        v11 = HelixTools().calculate_screw_axis(native)
        nativedegree = HelixTools().calc_angle(v11, 'z')

        diff = bestdegree - nativedegree
        tiltlist.append(diff)


### Build dataframe and save data
df = pd.DataFrame(columns=['PDB ID', 'RMSD', 'Depth', 'Tilt'])
df['PDB ID'] = proteinlist
df['RMSD'] = rmsdlist
df['Depth'] = depthlist
df['Tilt'] = tiltlist
df.to_csv('data.csv', sep=',')


