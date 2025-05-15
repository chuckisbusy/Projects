from pyrosetta import *
from modules import AddSpanlessMembraneMover
import os 
import re
import argparse
import numpy as np

np.set_printoptions(legacy='1.25')
init()
pymol = PyMOLMover()

show_in_pymol = True

###calculate the accuracy of the results####
def Intersection(result_embedded, start_embedded):
    TP = [value for value in result_embedded if value in start_embedded]
    FP = [value for value in result_embedded if value not in start_embedded]
    TN = [value for value in res_range if value not in result_embedded and value not in start_embedded]
    FN = [value for value in start_embedded if value not in result_embedded]

    TN = [str(x) for x in TN]
    
    # Remove any empty strings
    for i in (TP, FP, TN, FN):
        for j in i[:]:
            if not j.strip():
                i.remove(j)

    tp = len(TP)
    fp = len(FP)
    tn = len(TN)
    fn = len(FN)
    
    print("tp, fp, tn, fn")
    print(tp, fp, tn, fn)
    
    numerator = (tp * tn) - (fp * fn)
    denom = (tp+fp) * (tp+fn) * (tn+fp) * (tn+fn)
    
    if denom == 0:
        print("MCC cannot be calculated.")
        with open("results/{}_{}/txt/check_accuracy_AT.txt".format(protein_tag,multiple_tag), 'a+') as check_accuracy:
            check_accuracy.write("*" * 30)
            check_accuracy.write("\nThe results for {} {}\n".format(protein.split(".", 1)[0], y))
            check_accuracy.write("True positives: {}\n".format(TP))
            check_accuracy.write("False positives: {}\n".format(FP))
            check_accuracy.write("True negatives: {}\n".format(TN))
            check_accuracy.write("False negatives: {}\n".format(FN))
            check_accuracy.write("The MCC cannot be calulated.\n")
            check_accuracy.write("*" * 30)
            check_accuracy.close()
        sys.exit()
    elif denom != 0:
        pass
    
    MCC = (numerator) / (pow(denom,0.5))
    print(MCC)

    if len(start_embedded) and len(result_embedded) != 0:
        positive_perc = (len(TP) / len(start_embedded)) * 100
        false_positive_perc = (len(FP) / len(result_embedded)) * 100
        
        with open("results/{}_{}/txt/check_accuracy_AT.txt".format(protein_tag,multiple_tag), 'a+') as check_accuracy:
            check_accuracy.write("*" * 30)
            check_accuracy.write("\nThe results for {} {}\n".format(protein.split(".", 1)[0], y))
            check_accuracy.write("True positives: {}\n".format(TP))
            check_accuracy.write("False positives: {}\n".format(FP))
            check_accuracy.write("True negatives: {}\n".format(TN))
            check_accuracy.write("False negatives: {}\n".format(FN))
            check_accuracy.write("Positive percentage: {}%\n".format(positive_perc))
            check_accuracy.write("False positive percentage: {}%\n".format(false_positive_perc))
            check_accuracy.write("The MCC is: {}%\n".format(MCC))
            check_accuracy.write("*" * 30)

parser = argparse.ArgumentParser()
parser.add_argument('-protein',dest='protein',type=str,help="Name of input pdb")
parser.add_argument('-thickness', dest='thickness', type=float, default=15, help="Membrane thickness")
args = parser.parse_args()

#start with the OPM structure as the reference
protein = args.protein
protein_tag = protein.split(sep='_')[0]
multiple_tag = protein.split(sep='_')[1]
thickness = args.thickness

#delete old output files
if os.path.exists("results/{}_{}/txt/embedded_all_AT.txt".format(protein_tag,multiple_tag)):
    os.remove("results/{}_{}/txt/embedded_all_AT.txt".format(protein_tag,multiple_tag))
if os.path.exists("results/{}_{}/txt/embedded_opm_AT.txt".format(protein_tag,multiple_tag)):
    os.remove("results/{}_{}/txt/embedded_opm_AT.txt".format(protein_tag,multiple_tag))
if os.path.exists("results/{}_{}/txt/check_accuracy_AT.txt".format(protein_tag,multiple_tag)):
    os.remove("results/{}_{}/txt/check_accuracy_AT.txt".format(protein_tag,multiple_tag))

y = 35

#native pose
pose = pose_from_pdb("input_pdbs/{}_{}_renum.pdb".format(protein_tag,multiple_tag))
Pose = pyrosetta.rosetta.core.pose.Pose(pose)
        
#initiate the membrane object
mover = AddSpanlessMembraneMover()
mover.add_membrane_virtual(Pose)
mover.apply(Pose)
if show_in_pymol is True:
    pymol.apply(Pose)
        
#look for residues with atoms at least 1A deep inside the membrane
position = []
calculated_embedded = []
clean_duplicates = []
for i in range(1, Pose.size() - 1):
    j = Pose.residue(i).natoms()
    for atom in range(1, j+1):
        name = Pose.residue(i).atom_name(atom)
        position.append(Pose.residue(i).xyz(name).z)
        if Pose.residue(i).xyz(name).z > -thickness+1:
            calculated_embedded.append(i)
            clean_duplicates = []
            for k in range(0,len(calculated_embedded)):
                if k == 0:
                    clean_duplicates.append(calculated_embedded[k])
                elif k + 1 >= len(calculated_embedded):
                    break
                elif calculated_embedded[k+1] != calculated_embedded[k]:
                    clean_duplicates.append(calculated_embedded[k+1])
                else:
                    continue
with open("results/{}_{}/txt/embedded_opm_AT.txt".format(protein_tag,multiple_tag), 'a+') as record:
    record.write("The embedded residues are: {}\n".format(clean_duplicates))

#create a separate embedded data file for each scanned value for the result structures

if os.path.exists("results/{}_{}/txt/embedded_{}_{}_AT.txt".format(protein_tag,multiple_tag,protein_tag,multiple_tag)):
    os.remove("results/{}_{}/txt/embedded_{}_{}_AT.txt".format(protein_tag,multiple_tag,protein_tag,multiple_tag))

pose = pose_from_pdb("results/{}_{}/output_pdbs/{}_{}_renum_best_pose_overall.pdb".format(protein_tag,multiple_tag,protein_tag,multiple_tag))
Pose = pyrosetta.rosetta.core.pose.Pose(pose)

mover = AddSpanlessMembraneMover()
mover.add_membrane_virtual(Pose)
mover.apply(Pose)
mover.thickness = thickness
if show_in_pymol is True:
    pymol.apply(Pose)

#check whether the final geometry is in the positive or negative face
position = []
positive_face = ""
for i in range(1, Pose.size() - 1):
    j = Pose.residue(i).natoms()
    for atom in range(1, j+1):
        name = Pose.residue(i).atom_name(atom)
        position.append(Pose.residue(i).xyz(name).z)

if any(position[i] >= -thickness for i in range(0,len(position))):
    positive_face = True
elif all(position[i] <= thickness for i in range(0,len(position))):
    positive_face = False
else:
    if Pose.residue(1).xyz(1).z > 0:
        positive_face = True
    elif Pose.residue(1).xyz(1).z < 0:
        positive_face = False
print(positive_face)
clean_duplicates = []
calculated_embedded = []
if positive_face is True:
    for i in range(1, Pose.size() - 1):
        j = Pose.residue(i).natoms()
        for atom in range(1, j+1):
            name = Pose.residue(i).atom_name(atom)
            if Pose.residue(i).xyz(name).z <= thickness-1:
                calculated_embedded.append(i)
                clean_duplicates = []
                for k in range(0,len(calculated_embedded)):
                    if k == 0:
                        clean_duplicates.append(calculated_embedded[k])
                    elif k + 1 >= len(calculated_embedded):
                        break
                    elif calculated_embedded[k+1] != calculated_embedded[k]:
                        clean_duplicates.append(calculated_embedded[k+1])
                    else:
                        continue
    
    if len(calculated_embedded) != 0:
        with open("results/{}_{}/txt/embedded_{}_{}_AT.txt".format(protein_tag,multiple_tag,protein_tag,multiple_tag), 'a+') as record:
            record.write("The embedded residues are: {}\n".format(clean_duplicates))
        with open("results/{}_{}/txt/embedded_all_AT.txt".format(protein_tag,multiple_tag), 'a+') as record:
            record.write("The embedded residues for {}_{} {} are: {}\n".format(protein_tag, multiple_tag , y , clean_duplicates))

elif positive_face is False:
    for i in range(1, Pose.size() - 1):
        j = Pose.residue(i).natoms()
        for atom in range(1, j+1):
            name = Pose.residue(i).atom_name(atom)
            if Pose.residue(i).xyz(name).z >= -thickness+1:
                calculated_embedded.append(i)
                clean_duplicates = []
                for k in range(0,len(calculated_embedded)):
                    if k == 0:
                        clean_duplicates.append(calculated_embedded[k])
                    elif k + 1 >= len(calculated_embedded):
                        break
                    elif calculated_embedded[k+1] != calculated_embedded[k]:
                        clean_duplicates.append(calculated_embedded[k+1])
                    else:
                        continue
    if len(calculated_embedded) != 0:                
        with open("results/{}_{}/txt/embedded_{}_{}_AT.txt".format(protein_tag,multiple_tag,protein_tag,multiple_tag), 'a+') as record:
            record.write("The embedded residues are: {}\n".format(clean_duplicates))
        with open("results/{}_{}/txt/embedded_all_AT.txt".format(protein_tag,multiple_tag), 'a+') as record:
            record.write("The embedded residues for {}_{} {} are: {}\n".format(protein_tag, multiple_tag , y , clean_duplicates))
else:
    print("There's something wrong.\nಠ_ಠ")

with open("results/{}_{}/txt/embedded_opm_AT.txt".format(protein_tag,multiple_tag), 'r') as start_file:
    start_embedded = start_file.read()
    start_embedded = start_embedded.split()
with open("results/{}_{}/txt/embedded_{}_{}_AT.txt".format(protein_tag,multiple_tag,protein_tag,multiple_tag), 'r') as result_file:
    result_embedded = result_file.read()
    result_embedded = result_embedded.split()

for i in range(0,len(start_embedded)):
    start_embedded[i] = re.sub('[^A-Za-z0-9]+', '', start_embedded[i])
for i in range(0,len(result_embedded)):    
    result_embedded[i] = re.sub('[^A-Za-z0-9]+', '', result_embedded[i])

start_embedded = start_embedded[4:]
result_embedded = result_embedded[4:]

res_range = np.arange(1,Pose.size())
Intersection(result_embedded, start_embedded)
