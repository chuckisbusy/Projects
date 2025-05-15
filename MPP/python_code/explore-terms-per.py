#!/usr/bin/env python

import sys
import os
import numpy as np
import pandas as pd
import pyrosetta
from pyrosetta import *
from pyrosetta.teaching import *
from pyrosetta import rosetta

init()

def get_sasa(pose):
    '''Calculate the total and hydrophobic sasa'''
    atom_sasa = pyrosetta.rosetta.core.id.AtomID_Map_double_t()
    rsd_sasa = pyrosetta.rosetta.utility.vector1_double()
    rsd_hydrophobic_sasa = pyrosetta.rosetta.utility.vector1_double()
    rosetta.core.scoring.calc_per_res_hydrophobic_sasa(pose, rsd_sasa, rsd_hydrophobic_sasa, 1.4)
    tot_atom_sasa = rosetta.core.scoring.calc_per_atom_sasa(pose, atom_sasa, rsd_sasa, 1.4)

    return sum(rsd_sasa), sum(rsd_hydrophobic_sasa), tot_atom_sasa


#### Import Data
df=pd.read_csv("/home/cadeniran/ipn/data/selectd.csv", header=0, sep=',', engine='python')
df_per = df.drop(df[df['type_id'] == 1].index)

### Calculate TM
data_path = "/home/cadeniran/u1"
total_sasa, hydrophobic_sasa, fa_sols, atm_sasa = [], [], [], []
for x in df_per['pdbid']:
    dfile = x+"_A.pdb"
    pose = rosetta.core.pose.Pose()
    rosetta.core.import_pose.pose_from_file(pose, os.path.join(data_path, dfile))
    a,b,c = get_sasa(pose)
    total_sasa.append(a)
    hydrophobic_sasa.append(b)
    atm_sasa.append(c)
    scfxn = ScoreFunction()
    scfxn.set_weight(fa_sol, 1.0)
    fa_sols.append(scfxn.score(pose))

df_per["Total_SA"] = total_sasa
df_per["Hydrophobic_SA"] = hydrophobic_sasa
df_per["Total_Atom_SA"] = atm_sasa
df_per["fa_sol"] = fa_sols



os.chdir("/home/cadeniran/ipn/data")
df_per.to_csv('peripheral.csv', sep=',')

