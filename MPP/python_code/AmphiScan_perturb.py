from pyrosetta import *
from modules import AddSpanlessMembraneMover
from modules import HelixTools
from modules import HydrophobicMoment
import pyrosetta.rosetta.protocols.rigid as rigid_moves
import numpy as np
import os
import shutil
import sys
import timeit
import pandas as pd
import random as rand
import argparse
start = timeit.default_timer()

init()

pymol = PyMOLMover()

#Random Perturbation Movers
rotate = pyrosetta.rosetta.protocols.rigid.WholeBodyRotationMover
translate = pyrosetta.rosetta.protocols.rigid.WholeBodyTranslationMover
Vector = pyrosetta.rosetta.numeric.xyzVector_double_t

#Arguments for the code
parser = argparse.ArgumentParser()
parser.add_argument('-input_pdb', dest='input_pdb', type=str, help="Name of the input pdb.")
parser.add_argument('-pert_num', dest='pert_num', type=int, default=15, help="Number of perturbations performed at each depth. (default: 5)")
parser.add_argument('-z_increment', dest='z_increment', type=float, default=1, help="Increment of movement along the z-axis. default: 1 Angstrom")
args = parser.parse_args()

protein = args.input_pdb
protein_tag = protein.split(sep='_')[0]
multiple_tag = protein.split(sep='_')[1]

#create the folders to output the results
if not os.path.exists('results'):
    os.makedirs('results')
if not os.path.exists('results/{}_{}'.format(protein_tag,multiple_tag)):
    os.makedirs('results/{}_{}'.format(protein_tag,multiple_tag))
if not os.path.exists('results/{}_{}/txt'.format(protein_tag,multiple_tag)):
    os.makedirs('results/{}_{}/txt'.format(protein_tag,multiple_tag))
if not os.path.exists('results/{}_{}/csv'.format(protein_tag,multiple_tag)):
    os.makedirs('results/{}_{}/csv'.format(protein_tag,multiple_tag))
if not os.path.exists('results/{}_{}/output_pdbs'.format(protein_tag,multiple_tag)):
    os.makedirs('results/{}_{}/output_pdbs'.format(protein_tag,multiple_tag))

#load the pose based on the input value
Pose1 = pose_from_pdb('input_pdbs/{}_{}_renum.pdb'.format(protein_tag,multiple_tag))

thickness_from_file = 15
#scan parameters  
z_final = 25
starting_z = -26
z_increment = args.z_increment

#add the mp hbond option
pyrosetta.rosetta.core.scoring.hbonds.HBondOptions().mphbond()

#initiate the spanless membrane mover fom the modules file
fm = AddSpanlessMembraneMover()
fm.thickness = thickness_from_file
fm.membrane_core = thickness_from_file
fm.add_membrane_virtual(Pose1)
fm.apply(Pose1)

#Create a clone for the pose
rePose = Pose1.clone()

#Helix size
range_start = 1
range_end = rePose.size()-1

#calculate the center of mass of the protein
cmass = pyrosetta.rosetta.core.pose.center_of_mass(rePose, range_start, range_end)

#move the helix to the center of mass (0,0,0)
realign1 = translate(cmass.negated())
realign1.apply(rePose)

#move the helix center to the starting distance if one was set
if starting_z != 0:
    realign2 = translate(Vector(0,0,starting_z))
    realign2.apply(rePose)

#set the membrane scoring function
mem_sfxn = pyrosetta.rosetta.core.scoring.ScoreFunction()
mem_sfxn.add_weights_from_file("mpframework_smooth_fa_2012.wts")

#define the scan region
distance_from_membrane = 10
thickness = thickness_from_file
scan_start_region = -starting_z - thickness - distance_from_membrane
scan_stop_region = -starting_z + thickness + distance_from_membrane

#define Perturbation Movers
movemap = MoveMap()
movemap.set_bb(False)
movemap.set_chi(False)
movemap.set_jump(1,True)
kT = 1.0
nmoves = 1

#set starting values
best_z = 0
score_best = 999999
best_pose = Pose()
local_best = Pose()
rePose2 = Pose()
rePose3 = Pose()

#start the scan along the z-axis
zRange = np.arange(scan_start_region, scan_stop_region+1, z_increment)

numPert = args.pert_num

rePose1 = rePose.clone()

for dist in zRange:
    if os.path.exists("results/{}_{}/txt/{}_scores_{}.txt".format(protein_tag,multiple_tag, protein.split(".",1)[0], dist)):
        os.remove("results/{}_{}/txt/{}_scores_{}.txt".format(protein_tag,multiple_tag, protein.split(".",1)[0], dist))
    if os.path.exists("results/{}_{}/csv/{}_socres_{}.csv".format(protein_tag,multiple_tag, protein.split(".",1)[0], dist)):
        os.remove("results/{}_{}/csv/{}_scores_{}.csv".format(protein_tag,multiple_tag, protein.split(".",1)[0], dist))

    rePose2 = rePose1.clone()

    #Move down z-axis
    trans_z = translate(Vector(0,0,dist))
    trans_z.apply(rePose2)

    #create the spin mover
    spin_mover = pyrosetta.rosetta.protocols.rigid.RigidBodySpinMover(1)
    spin_mover.rot_center(pyrosetta.rosetta.core.pose.center_of_mass(rePose2, range_start, range_end))
    spin_mover.spin_mag(60.0)
    spin_mover.spin_axis(Vector(0,1,0))
    spin_mover.apply(rePose2)
    
    score_local_best = 99999
    local_pose = rePose2.clone()

    #use a 3-step procedure to zero-in to the best pose. The rotation perturbation amount goes down with each step by half 20/10/5
    for i in range(1,numPert+1):
        rePose3 = rePose2.clone()

        # Random move with RandomizeMover
        randomize1 = pyrosetta.rosetta.protocols.rigid.RigidBodyRandomizeMover(rePose3, 1, rigid_moves.partner_upstream)
        randomize1.apply(rePose3)

        # Monte Carlo Spin movement
        mhm = pyrosetta.rosetta.protocols.canonical_sampling.MetropolisHastingsMover()
        spin_mover.rot_center(pyrosetta.rosetta.core.pose.center_of_mass(rePose3, range_start, range_end))
        spin_mover.spin_mag(8.0)
        spin_mover.initialize_simulation(rePose3, mhm, 8000)
        spin_mover.apply(rePose3)

        # Minimize
        min_mover = rosetta.protocols.minimization_packing.MinMover()
        min_mover.movemap(movemap)
        min_mover.tolerance(1e-6)
        min_mover.apply(rePose3)

        score_new = mem_sfxn(rePose3)

        if score_new < score_local_best:
            score_local_best = score_new
            local_pose = rePose3.clone()

        if score_new < score_best:
            best_z = dist
            best_pose = rePose3.clone()
            best_iteration = i
            score_best = score_local_best

        if i == numPert:
            #repack the best pose
            task_pack = pyrosetta.rosetta.core.pack.task.TaskFactory.create_packer_task(local_pose)
            task_pack.restrict_to_repacking()
            membrane_pack = pyrosetta.rosetta.protocols.minimization_packing.PackRotamersMover(mem_sfxn, task_pack)
            membrane_pack.apply(local_pose)

            score_local_best = mem_sfxn(local_pose)

            if score_local_best < score_best:
                best_z = dist
                best_pose = local_pose.clone()
                score_best = score_local_best

            #record the best score of all perturbations at all z-values to txt and csv
            with open("results/{}_{}/txt/{}_best_scores.txt".format(protein_tag,multiple_tag, protein.split(".",1)[0]), "a+") as best:
                best.write("The best score of {} belongs to {} Angstroms.\n".format(score_local_best, dist))

            with open("results/{}_{}/csv/{}_best_scores.csv".format(protein_tag,multiple_tag, protein.split(".",1)[0]), "a+") as best:
                best.write("{}, {}\n".format(dist, score_local_best))

            #Dump best pose at current depth
            local_pose.dump_pdb("results/{}_{}/output_pdbs/{}_{}_pose_{}.pdb".format(protein_tag,multiple_tag, protein.split(".",1)[0], "best",dist))

    if dist > scan_stop_region:
        best_pose.dump_pdb("results/{}_{}/output_pdbs/{}_{}_pose_{}.pdb".format(protein_tag,multiple_tag, protein.split(".",1)[0], "best","overall"))

        #record the overall best score and z-value
        with open("results/{}_{}/txt/{}_best_scores.txt".format(protein_tag,multiple_tag, protein.split(".",1)[0]), "a+") as best:
            best.write("The best depth {} has score {}\n".format(best_z, score_best))
            best.write("Time elapsed for completion: {} seconds.\n".format(timeit.default_timer() - start))
