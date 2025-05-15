from pyrosetta import * 
from modules import AddSpanlessMembraneMover
from modules import HelixTools
from modules import HydrophobicMoment
from pyrosetta import PyMOLMover
import argparse
import numpy as np
import pandas as pd
import string
import os

init()
pymol = PyMOLMover()

parser = argparse.ArgumentParser()
parser.add_argument('-protein',dest='protein', type=str, help="Name of the protein file")
args = parser.parse_args()

protein = args.protein
protein_tag = protein.split(sep='_')[0]
multiple_tag = protein.split(sep='_')[1]

native = pose_from_pdb('input_pdbs/{}_{}_renum.pdb'.format(protein_tag,multiple_tag))
best = pose_from_pdb('results/{}_{}/output_pdbs/{}_{}_renum_best_pose_overall.pdb'.format(protein_tag,multiple_tag,protein_tag,multiple_tag))
pymol.apply(best)

#def AH_CA_rmsd(native, best):
number_of_residues = native.size()

#initiate the spanless membrane mover
fm = AddSpanlessMembraneMover()
fm.add_membrane_virtual(native)
fm.apply(native)
pymol.apply(native)

#calculate the center of mass for the best and native helices
cmass_best = pyrosetta.rosetta.core.pose.center_of_mass(best, 1, best.size()-1)
cmass_native = pyrosetta.rosetta.core.pose.center_of_mass(native, 1, native.size()-1)

print("The best cmass is: {}".format(cmass_best))
print("The native cmass is: {}".format(cmass_native))

move_xy = pyrosetta.rosetta.numeric.xyzVector_double_t(-cmass_native[0], -cmass_native[1],0)

#run the translation
shifted = native.clone()
copy_best = best.clone()
translation_mover = pyrosetta.rosetta.protocols.rigid.WholeBodyTranslationMover(move_xy)
translation_mover.apply(shifted)
pymol.apply(shifted)

print(f"shifted COM: {pyrosetta.rosetta.core.pose.center_of_mass(shifted, 1, shifted.size()-1)}")

#calculate the helix parameters
ht = HelixTools()
helix_normal = ht.calculate_screw_axis(shifted)
native_angle_with_x = ht.calc_angle(helix_normal,'x')
native_angle_with_y = ht.calc_angle(helix_normal,'y')
native_angle_with_z = ht.calc_angle(helix_normal,'z')
print(f"shifted: {native_angle_with_x}, {native_angle_with_y}, {native_angle_with_z}")

helix_best_normal = ht.calculate_screw_axis(best)
best_angle_with_x = ht.calc_angle(helix_best_normal,'x')
best_angle_with_y = ht.calc_angle(helix_best_normal,'y')
best_angle_with_z = ht.calc_angle(helix_best_normal,'z')
print(f"best: {best_angle_with_x}, {best_angle_with_y}, {best_angle_with_z}") 

x_diff = best_angle_with_x - native_angle_with_x
print(f"x diff: {x_diff}")

if cmass_best[2] > 0:
    if x_diff < 0:
        align_x = pyrosetta.rosetta.protocols.rigid.WholeBodyRotationMover(pyrosetta.rosetta.numeric.xyzVector_double_t(1,0,1), 
                                                                            pyrosetta.rosetta.core.pose.center_of_mass(shifted, 1, shifted.size()-1), 
                                                                            -(x_diff))
        align_x.apply(shifted)
    else:
        align_x = pyrosetta.rosetta.protocols.rigid.WholeBodyRotationMover(pyrosetta.rosetta.numeric.xyzVector_double_t(1,0,1), 
                                                                            pyrosetta.rosetta.core.pose.center_of_mass(shifted, 1, shifted.size()-1), 
                                                                            (x_diff))
        align_x.apply(shifted)
else:
    if x_diff > 0:
        align_x = pyrosetta.rosetta.protocols.rigid.WholeBodyRotationMover(pyrosetta.rosetta.numeric.xyzVector_double_t(1,0,1), 
                                                                            pyrosetta.rosetta.core.pose.center_of_mass(shifted, 1, shifted.size()-1), 
                                                                            -(x_diff))
        align_x.apply(shifted)
    else:
        align_x = pyrosetta.rosetta.protocols.rigid.WholeBodyRotationMover(pyrosetta.rosetta.numeric.xyzVector_double_t(1,0,1), 
                                                                            pyrosetta.rosetta.core.pose.center_of_mass(shifted, 1, shifted.size()-1), 
                                                                            (x_diff)) ### 4nwz
        align_x.apply(shifted)
pymol.apply(shifted)

helix_normal = ht.calculate_screw_axis(shifted)
native_angle_with_x = ht.calc_angle(helix_normal,'x')
native_angle_with_y = ht.calc_angle(helix_normal,'y')
native_angle_with_z = ht.calc_angle(helix_normal,'z')
print(f"shifted: {native_angle_with_x}, {native_angle_with_y}, {native_angle_with_z}")    

y_diff = best_angle_with_y - native_angle_with_y
print(f"y diff: {y_diff}")

if cmass_best[2] > 0:
    if y_diff < 0:
        align_y = pyrosetta.rosetta.protocols.rigid.WholeBodyRotationMover(pyrosetta.rosetta.numeric.xyzVector_double_t(1,0,1), 
                                                                            pyrosetta.rosetta.core.pose.center_of_mass(shifted, 1, shifted.size()-1), 
                                                                            (y_diff))
        align_y.apply(shifted)
    else:
        align_y = pyrosetta.rosetta.protocols.rigid.WholeBodyRotationMover(pyrosetta.rosetta.numeric.xyzVector_double_t(1,0,1), 
                                                                            pyrosetta.rosetta.core.pose.center_of_mass(shifted, 1, shifted.size()-1), 
                                                                            -(y_diff))
        align_y.apply(shifted)
else:
    if y_diff > 0:
        align_y = pyrosetta.rosetta.protocols.rigid.WholeBodyRotationMover(pyrosetta.rosetta.numeric.xyzVector_double_t(1,0,1), 
                                                                            pyrosetta.rosetta.core.pose.center_of_mass(shifted, 1, shifted.size()-1), 
                                                                            (y_diff)) ### 4nwz
        align_y.apply(shifted)
    else:
        align_y = pyrosetta.rosetta.protocols.rigid.WholeBodyRotationMover(pyrosetta.rosetta.numeric.xyzVector_double_t(1,0,1), 
                                                                            pyrosetta.rosetta.core.pose.center_of_mass(shifted, 1, shifted.size()-1), 
                                                                            -(y_diff))
        align_y.apply(shifted)
pymol.apply(shifted)

helix_native_normal = ht.calculate_screw_axis(shifted)
native_angle_with_x = ht.calc_angle(helix_native_normal,'x')
native_angle_with_y = ht.calc_angle(helix_native_normal,'y')
native_angle_with_z = ht.calc_angle(helix_native_normal,'z')


x_diff = best_angle_with_x - native_angle_with_x
y_diff = best_angle_with_y - native_angle_with_y
z_diff = best_angle_with_z - native_angle_with_z
print(f"xyz diff: {x_diff}, {y_diff}, {z_diff}")
shifted.dump_pdb("shifted_native.pdb")
pymol.apply(shifted)

ca_shifted = []
ca_best = []

#calculate the CA positions for the shifted original and the best structures
for i in range(1,shifted.size()):
    ca_shifted.append(np.array(shifted.residue(i).xyz('CA')))
for j in range(1, best.size()):
    ca_best.append(np.array(copy_best.residue(j).xyz('CA')))

shifted_df = pd.DataFrame(ca_shifted)
best_df = pd.DataFrame(ca_best)

total = 0
for k in range(number_of_residues):
    if all(item < 0 for item in shifted_df[2]) is all(item < 0 for item in best_df[2]):
        total = total + np.square((ca_shifted[k][0] - ca_best[k][0])) + np.square((ca_shifted[k][1] - ca_best[k][1])) + np.square((ca_shifted[k][2] - ca_best[k][2]))
    else:
        total = total + np.square((abs(ca_shifted[k][0]) - abs(ca_best[k][0]))) + np.square((abs(ca_shifted[k][1]) - abs(ca_best[k][1]))) + np.square((abs(ca_shifted[k][2]) - abs(ca_best[k][2])))

average = total/len(ca_best)
squared_average = np.sqrt(average)

print(f"Avg: {squared_average}")
print(squared_average)

with open('results/{}_{}/txt/rmsd.txt'.format(protein_tag,multiple_tag), 'w') as rmsd_file:
    rmsd_file.write("The RMSD is: {}".format(squared_average))

with open('corrdata.csv', 'a') as xyz_rmsd:
    xyz_rmsd.write(f"{protein}, {x_diff}, {y_diff}, {z_diff}, {squared_average} \n")
    xyz_rmsd.close()
