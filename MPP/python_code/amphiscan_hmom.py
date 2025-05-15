from pyrosetta import *
from modules import AddSpanlessMembraneMover
from modules import HelixTools
from modules import HydrophobicMoment
from modules import HydrophCalc
import modules
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
parser.add_argument('-pert_num', dest='pert_num', type=int, default=5, help="Number of perturbations performed at each depth. (default: 5)")
parser.add_argument('-z_increment', dest='z_increment', type=int, default=1, help="Increment of movement along the z-axis. default: 1 Angstrom")
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
cmass = HydrophCalc(Pose1)

thickness_from_file = 15
#scan parameters  
z_final = 25
starting_z = -26
z_increment = 1



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
#cmass = pyrosetta.rosetta.core.pose.center_of_mass(rePose, range_start, range_end)

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
movemap.set_bb(True)
kT = 1.0
nmoves = 1

#set starting values
best_z = 0
score_best = 999999
best_pose = Pose()
rePose2 = Pose()
rePose3 = Pose()

#start the scan along the z-axis
zRange = np.arange(scan_start_region, scan_stop_region+z_increment, z_increment)

numPert = args.pert_num

#rePose2 = rePose.clone()

#"""
for dist in zRange:

    if os.path.exists("results/{}_{}/txt/{}_scores_{}.txt".format(protein_tag,multiple_tag, protein.split(".",1)[0], dist)):
        os.remove("results/{}_{}/txt/{}_scores_{}.txt".format(protein_tag,multiple_tag, protein.split(".",1)[0], dist))
    if os.path.exists("results/{}_{}/csv/{}_socres_{}.csv".format(protein_tag,multiple_tag, protein.split(".",1)[0], dist)):
        os.remove("results/{}_{}/csv/{}_scores_{}.csv".format(protein_tag,multiple_tag, protein.split(".",1)[0], dist))

    rePose2 = rePose.clone()

    #Move down z-axis
    trans_z = translate(Vector(0,0,dist))
    trans_z.apply(rePose2)

    rePose3 = rePose2.clone()

    score_local_best = 99999

    #use a 3-step procedure to zero-in to the best pose. The rotation perturbation amount goes down with each step by half 20/10/5
    for i in range(1,numPert+1):
        rePose3 = rePose2.clone()
        rotator = pyrosetta.rosetta.protocols.rigid.WholeBodyRotationMover(Vector(rand.randint(-1,1), rand.randint(-1,1), rand.randint(-1,1)),
                pyrosetta.rosetta.core.pose.center_of_mass(rePose3, range_start, range_end), 
                20)
        rotator.reinitialize_for_new_input()
        rotator.apply(rePose3)

        #Score Pose
        try:
            score_new = mem_sfxn(rePose3)
        except RuntimeError:
            print("ERROR:NAN occurred in H-bonding calculations!")

        if score_new < score_local_best:
            score_local_best = score_new
            rePose2 = rePose3.clone()

        if score_new < score_best:
            best_z = dist
            best_pose = rePose3.clone()
            best_iteration = i
            score_best = score_local_best

        #perturbations with a 10 degree rotation
        rotator = pyrosetta.rosetta.protocols.rigid.WholeBodyRotationMover(Vector(rand.randint(-1,1), rand.randint(-1,1), rand.randint(-1,1)),
                pyrosetta.rosetta.core.pose.center_of_mass(rePose3, range_start, range_end), 
                10)
        rotator.reinitialize_for_new_input()
        rotator.apply(rePose3)

        #Score Pose
        try:
            score_new = mem_sfxn(rePose3)
        except RuntimeError:
            print("ERROR:NAN occurred in H-bonding calculations!")

        if score_new < score_local_best:
            score_local_best = score_new
            rePose2 = rePose3.clone()

        if score_new < score_best:
            best_z = dist
            best_pose = rePose3.clone()
            best_iteration = i
            score_best = score_local_best

        #perturbations with a 5 degree rotation
        rotator = pyrosetta.rosetta.protocols.rigid.WholeBodyRotationMover(Vector(rand.randint(-1,1), rand.randint(-1,1), rand.randint(-1,1)),
                pyrosetta.rosetta.core.pose.center_of_mass(rePose3, range_start, range_end), 
                5)
        rotator.reinitialize_for_new_input()
        rotator.apply(rePose3)

        #Score Pose
        try:
            score_new = mem_sfxn(rePose3)
        except RuntimeError:
            print("ERROR:NAN occurred in H-bonding calculations!")

        if score_new < score_local_best:
            score_local_best = score_new
            rePose2 = rePose3.clone()

        if score_new < score_best:
            best_z = dist
            best_pose = rePose3.clone()
            best_iteration = i
            score_best = score_local_best

        #perturbations with a 1 degree rotation
        rotator = pyrosetta.rosetta.protocols.rigid.WholeBodyRotationMover(Vector(rand.randint(-1,1), rand.randint(-1,1), rand.randint(-1,1)),
                pyrosetta.rosetta.core.pose.center_of_mass(rePose3, range_start, range_end), 
                1)
        rotator.reinitialize_for_new_input()
        rotator.apply(rePose3)

        #Score Pose
        try:
            score_new = mem_sfxn(rePose3)
        except RuntimeError:
            print("ERROR:NAN occurred in H-bonding calculations!")

        if score_new < score_local_best:
            score_local_best = score_new
            rePose2 = rePose3.clone()

        if score_new < score_best:
            best_z = dist
            best_pose = rePose3.clone()
            best_iteration = i
            score_best = score_local_best

        #record all three perturbation's score to txt for this z-value
        with open("results/{}_{}/txt/{}_scores_{}.txt".format(protein_tag,multiple_tag, protein.split(".",1)[0], dist), "a+") as all_scores:
            all_scores.write("The score of {} belongs to {} Angstroms.\n".format(score_new, dist))
        #record all three perturbation's score to csv for this z-value 
        with open("results/{}_{}/csv/{}_scores_{}.csv".format(protein_tag,multiple_tag, protein.split(".",1)[0], dist), "a+") as all_scores_csv:
                all_scores_csv.write("{}, {}\n".format(dist, score_new))

        if i == numPert:
             #record the best score of all perturbations made to txt for this z-value
            with open("results/{}_{}/txt/{}_scores_{}.txt".format(protein_tag,multiple_tag, protein.split(".",1)[0], dist), "a+") as all_scores:
               all_scores.write("The best score of {} belongs to {} Angstroms.\n".format(score_local_best, dist))

            #record the best score of all perturbations at all z-values to txt and csv
            with open("results/{}_{}/txt/{}_best_scores.txt".format(protein_tag,multiple_tag, protein.split(".",1)[0]), "a+") as best:
                best.write("The best score of {} belongs to {} Angstroms.\n".format(score_local_best, dist))

            with open("results/{}_{}/csv/{}_best_scores.csv".format(protein_tag,multiple_tag, protein.split(".",1)[0]), "a+") as best:
                best.write("{}, {}\n".format(dist, score_local_best))

            #Dump best pose at current depth
            rePose2.dump_pdb("results/{}_{}/output_pdbs/{}_{}_pose_{}.pdb".format(protein_tag,multiple_tag, protein.split(".",1)[0], "best",dist))

    if dist >= scan_stop_region:
        best_pose.dump_pdb("results/{}_{}/output_pdbs/{}_{}_pose_{}.pdb".format(protein_tag,multiple_tag, protein.split(".",1)[0], "best","overall"))
        print("#" * 30)
        print("DONE")
        print("#" * 30)
        #record the overall best score and z-value
        with open("results/{}_{}/txt/{}_best_scores.txt".format(protein_tag,multiple_tag, protein.split(".",1)[0]), "a+") as best:
            best.write("The best depth {} has score {}\n".format(best_z, score_best))
            best.write("Time elapsed for completion: {} seconds.\n".format(timeit.default_timer() - start))

#print("Completed in: {} seconds.\n".format(timeit.default_timer() - start))
