from pyrosetta import * 
from modules import AddSpanlessMembraneMover
from modules import HelixTools
from modules import HydrophobicMoment
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os

init()
pymol = PyMOLMover()


def get_tilt_angle(best_pose):
    """ Calculates the tilt angle of a pose by using existing methods in
        the HelixTools class. Must be run on a membrane-initiated pose. """
    ht = HelixTools()
    v1 = ht.calculate_screw_axis(best_pose)
    degree = ht.calc_angle(v1, 'z')
    return degree
    

def count_charged_residues(pose):
    counter = 0
    for i in range(1, len(pose.sequence())):
        if pose.residue(i).is_charged() == True:
            counter += 1
        else:
            continue
    return counter / pose.size()


def count_positively_charged_residues(pose):
    counter = 0
    for i in range(1, len(pose.sequence())):
        if (pose.residue(i).is_charged() == True) and (pose.sequence(i,i) in ['K', 'R']):
            counter += 1
        else:
            continue
    return counter / pose.size()


def count_negatively_charged_residues(pose):
    counter = 0
    for i in range(1, len(pose.sequence())):
        if (pose.residue(i).is_charged() == True) and (pose.sequence(i,i) in ['D', 'E']):
            counter += 1
        else:
            continue
    return counter / pose.size()


def count_polar_residues(pose):
    counter = 0
    for i in range(1, len(pose.sequence())):
        if pose.residue(i).is_polar() == True:
            counter += 1
        else:
            continue
    return counter / pose.size()


def count_apolar_residues(pose):
    counter = 0
    for i in range(1, len(pose.sequence())):
        if pose.residue(i).is_apolar() == True:
            counter += 1
        else:
            continue
    return counter / pose.size()


def count_aliphatic_residues(pose):
    counter = 0
    for i in range(1, len(pose.sequence())):
        if pose.sequence(i,i) in ['G', 'A', 'V', 'L', 'I']:
            counter += 1
        else:
            continue
    return counter / pose.size()


def count_aromatic_residues(pose):
    counter = 0
    for i in range(1, len(pose.sequence())):
        if pose.residue(i).is_aromatic() == True:
            counter += 1
        else:
            continue
    return counter / pose.size()


def calculate_embedded_residues(protein_name, thickness=15):
    """ This function takes the name of the protein and calculates the embedded residues.
        calculated as any residue that has at least one atom within 1A of the
        membrane. The percentage of positive and negative results are calculated and
        printed to a separate file. """
    #delete old output files
    if os.path.exists("results/{}/txt/embedded_all.txt".format(protein_name)):
        os.remove("results/{}/txt/embedded_all.txt".format(protein_name))
    if os.path.exists("results/{}/txt/embedded_opm.txt".format(protein_name)):
        os.remove("results/{}/txt/embedded_opm.txt".format(protein_name))
    if os.path.exists("results/{}/txt/check_accuracy.txt".format(protein_name)):
        os.remove("results/{}/txt/check_accuracy.txt".format(protein_name))

    native_pose = pose_from_pdb("input_pdbs/{}.pdb".format(protein_name))
    
    #initiate the membrane object
    membrane_mover = AddSpanlessMembraneMover()
    membrane_mover.add_membrane_virtual(native_pose)
    membrane_mover.thickness = thickness
    membrane_mover.apply(native_pose)

    #check whether the native geometry is in the positive or negative face
    native_position = []
    native_positive_face = ""
    for i in range(1, native_pose.size() - 1):
        j = native_pose.residue(i).natoms()
        for atom in range(1, j+1):
            name = native_pose.residue(i).atom_name(atom)
            native_position.append(native_pose.residue(i).xyz(name).z)
    
    if all(native_position[i] > -thickness for i in range(0,len(native_position))):
        native_positive_face = True
    elif all(native_position[i] < thickness for i in range(0,len(native_position))):
        native_positive_face = False
    
    #look for native residues with atoms at least 1A deep inside the membrane
    calculated_embedded = []
    clean_duplicates = []
    if native_positive_face is True:
        for i in range(1, native_pose.size() - 1):
            j = native_pose.residue(i).natoms()
            for atom in range(1, j+1):
                name = native_pose.residue(i).atom_name(atom)
                native_position.append(native_pose.residue(i).xyz(name).z)
                if native_pose.residue(i).xyz(name).z < thickness-1:
                    calculated_embedded.append(i)
                    clean_duplicates = []
                    for k in range(0,len(calculated_embedded)):
                        if k == 0:
                            clean_duplicates.append(calculated_embedded[k])
                        elif k + 1 >= len(calculated_embedded):
                            break
                        elif calculated_embedded[k+1] is not calculated_embedded[k]:
                            clean_duplicates.append(calculated_embedded[k+1])
                        else:
                            continue
                        
    elif native_positive_face is False:
        for i in range(1, native_pose.size() - 1):
            j = native_pose.residue(i).natoms()
            for atom in range(1, j+1):
                name = native_pose.residue(i).atom_name(atom)
                native_position.append(native_pose.residue(i).xyz(name).z)
                if native_pose.residue(i).xyz(name).z > -thickness+1:
                    calculated_embedded.append(i)
                    clean_duplicates = []
                    for k in range(0,len(calculated_embedded)):
                        if k == 0:
                            clean_duplicates.append(calculated_embedded[k])
                        elif k + 1 >= len(calculated_embedded):
                            break
                        elif calculated_embedded[k+1] is not calculated_embedded[k]:
                            clean_duplicates.append(calculated_embedded[k+1])
                        else:
                            continue
                    
    start_embedded = clean_duplicates            
    with open('results/{}/txt/embedded_opm.txt'.format(protein_name), 'a+') as record:
        record.write("The embedded residues are: {}\n".format(clean_duplicates))
        
    #create a separate embedded data file for each scanned value for the result structures
    if os.path.exists("results/{}/txt/embedded_all.txt".format(protein_name)):
        os.remove("results/{}/txt/embedded_all.txt".format(protein_name))

    best_pose = pose_from_pdb("results/{}/output_pdbs/best_{}_0.pdb".format(protein_name, protein_name))
    membrane_mover.apply(best_pose)

    #check whether the final geometry is in the positive or negative face
    position = []
    positive_face = ""
    for i in range(1, best_pose.size() - 1):
        j = best_pose.residue(i).natoms()
        for atom in range(1, j+1):
            name = best_pose.residue(i).atom_name(atom)
            position.append(best_pose.residue(i).xyz(name).z)
    
    if all(position[i] > -thickness for i in range(0,len(position))):
        positive_face = True
    elif all(position[i] < thickness for i in range(0,len(position))):
        positive_face = False
    
    calculated_embedded = []
    result_embedded = []
    if positive_face is True:
        for i in range(1, best_pose.size() - 1):
            j = best_pose.residue(i).natoms()
            for atom in range(1, j+1):
                name = best_pose.residue(i).atom_name(atom)
                if best_pose.residue(i).xyz(name).z < thickness-1:
                    calculated_embedded.append(i)
                    clean_duplicates = []
                    for k in range(0,len(calculated_embedded)):
                        if k == 0:
                            clean_duplicates.append(calculated_embedded[k])
                        elif k + 1 >= len(calculated_embedded):
                            break
                        elif calculated_embedded[k+1] is not calculated_embedded[k]:
                            clean_duplicates.append(calculated_embedded[k+1])
                        else:
                            continue
       
        if len(calculated_embedded) != 0:
            result_embedded = clean_duplicates                
            with open('results/{}/txt/embedded_all.txt'.format(protein_name), 'a+') as record:
                record.write("The embedded residues are: {}\n".format(clean_duplicates))
          
    elif positive_face is False:
        for i in range(1, best_pose.size() - 1):
            j = best_pose.residue(i).natoms()
            for atom in range(1, j+1):
                name = best_pose.residue(i).atom_name(atom)
                if best_pose.residue(i).xyz(name).z > -thickness+1:
                    calculated_embedded.append(i)
                    clean_duplicates = []
                    for k in range(0,len(calculated_embedded)):
                        if k == 0:
                            clean_duplicates.append(calculated_embedded[k])
                        elif k + 1 >= len(calculated_embedded):
                            break
                        elif calculated_embedded[k+1] is not calculated_embedded[k]:
                            clean_duplicates.append(calculated_embedded[k+1])
                        else:
                            continue
        if len(calculated_embedded) != 0:
            result_embedded = clean_duplicates                
            with open('results/{}/txt/embedded_all.txt'.format(protein_name), 'a+') as record:
                record.write("The embedded residues are: {}\n".format(clean_duplicates))
           
    else:
        print("There's something wrong.\nಠ_ಠ")
        sys.exit()
    
    res_range = np.arange(1,native_pose.size())
    true_positives = [value for value in result_embedded if value in start_embedded]
    false_positives = [value for value in result_embedded if value not in start_embedded]
    true_negatives = [value for value in res_range if value not in result_embedded and value not in start_embedded]
    false_negatives = [value for value in start_embedded if value not in result_embedded]
    
    #calculate MCC 
    tp = len(true_positives)
    fp = len(false_positives)
    tn = len(true_negatives)
    fn = len(false_negatives)
    mcc = (tp * tn - fp * fn) / np.sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
    
    if len(start_embedded) and len(result_embedded) is not 0:
        coverage = (len(true_positives) / len(start_embedded)) * 100
        positive_perc = (len(true_positives) / len(result_embedded)) * 100
        false_positive_perc = (len(false_positives) / len(result_embedded)) * 100
        accuracy = 100 * (len(true_positives) + len(true_negatives)) / len(res_range)
        
        with open('results/{}/txt/check_accuracy.txt'.format(protein_name), 'a+') as check_accuracy:
            check_accuracy.write("*" * 30)
            check_accuracy.write("\nThe results for {}\n".format(protein_name))
            check_accuracy.write("True positive residues: {}\n".format(true_positives))
            check_accuracy.write("False positive residues: {}\n".format(false_positives))
            check_accuracy.write("True negative residues: {}\n".format(true_negatives))
            check_accuracy.write("False negative residues: {}\n".format(false_negatives))
            check_accuracy.write("Query coverage: {}\n".format(round(coverage, 0)))
            check_accuracy.write("True positive percentage: {} %\n".format(round(positive_perc, 0)))
            check_accuracy.write("False positive percentage: {} %\n".format(round(false_positive_perc, 0)))
            check_accuracy.write("The accuracy is: {} %\n".format(round(accuracy,0)))
            check_accuracy.write("The MCC is: {}\n".format(mcc))
            check_accuracy.write("*" * 30)
    
    elif len(start_embedded) and len(result_embedded) is 0:
        with open('results/{}/txt/check_accuracy.txt'.format(protein_name), 'a+') as check_accuracy:
            check_accuracy.write("*" * 30)
            check_accuracy.write("\nThe results for {}\n".format(protein_name))
            check_accuracy.write("True positives: []\n")
            check_accuracy.write("False positives: []\n")
            check_accuracy.write("Missing residues: []\n")
            check_accuracy.write("Query coverage: 0 %\n")
            check_accuracy.write("True positive percentage: 100 %\n")
            check_accuracy.write("False positive percentage: 0 %\n")
            check_accuracy.write("*" * 30)
            
    else:
        with open('results/{}/txt/check_accuracy.txt'.format(protein_name), 'a+') as check_accuracy:
            check_accuracy.write("*" * 30)
            check_accuracy.write("\nThe results for {}\n".format(protein_name))
            check_accuracy.write("True positive percentage: NA\n")
            check_accuracy.write("Query coverage: NA\n")
            check_accuracy.write("The embedded file is emtpy!\n")
            check_accuracy.write("*" * 30)         


#getter functions for the functions already available in Rosetta
def get_size(native_pose):
    """ Gets the size of a pose using Rosetta size function. Should be used 
        with the native pose rather than the membrane pose. """
    pose_size = native_pose.size()
    return pose_size
    

def get_residue_frequency(native_pose):
    """ Calculates the frequency of each residue within a given pose. """
    sequence = native_pose.sequence()
    aa_list = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P',
               'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
    frequency = {}
    for aa in aa_list:
        frequency[aa] = round(sequence.count(aa)/len(sequence), 2)
    return frequency


def get_charge(native_pose):
    """ Calculates the charge of a pose based on its sequence using default
        Rosetta protonation states. Should be used with the native pose 
        rather than the membrane pose."""
    aa_charge_table = {'D': -1, 'E': -1, 'R': 1, 'K': 1, 'H': 0, 'A': 0,
                       'I': 0, 'L': 0, 'M': 0, 'V': 0, 'F': 0, 'W': 0,
                       'Y': 0, 'N': 0, 'Q': 0, 'C': 0, 'S': 0, 'T': 0,
                       'G': 0, 'P': 0}
    sequence = native_pose.sequence()
    charge  = 0
    for residue in sequence:
        charge = charge + aa_charge_table[residue]
    return charge

    
def get_number_of_prolines(native_pose):
    """ Calculates the number of prolines. Should be used with the native pose 
        rather than the membrane pose. """
    prolines = 0
    sequence = native_pose.sequence()
    for residue in sequence:
        if residue == 'P':
            prolines += 1
    return prolines


def sns_heatmap(protein, y_angle):
    """ Creates a heatmap from the helix rotation angle, membrane depth,
        and score belonging to a particular y-angle value. Bar positions are
        set based on 0.5 A z step size. """
    df = pd.read_csv('csv/{}_{}_xz_scores.csv'.format(protein, y_angle), names=['X', 'Y', 'Z'])
    df.iloc[:,0] = df.iloc[:,0] - 100
    
    data_pivoted = df.pivot_table(index='X', columns='Y', values='Z')
    
    ax = sns.heatmap(data_pivoted, cmap='jet')
    #membrane positive
    ax.hlines(y=20.5, xmin=0, xmax=180, linewidth=2, color='black', linestyle='solid')
    ax.hlines(y=80.5, xmin=0, xmax=180, linewidth=2, color='black', linestyle='solid')
    #membrane negative
    ax.hlines(y=26.5, xmin=0, xmax=180, linewidth=1, color='blue', linestyle='dashed')
    ax.hlines(y=74.5, xmin=0, xmax=180, linewidth=1, color='blue', linestyle='dashed')
    
    plt.xlabel('x-angle (degrees)')
    plt.ylabel('Distance from membrane (A)')
    plt.savefig('{}_{}.png'.format(protein, y_angle), dpi=800)


def analyze_y():
    """ Plots the y-angles versus the best Rosetta score from each structure.
        It must be run in the txt folder. """

    """ Can merge the rmsd and heatmap parts under a single class that stores the 
        RMSD information and passes it to seaborn for plotting """
    y_values = []
    scores = []
    
    os.system("grep best txt/*.txt | awk '{print $1, $5}' > txt/y_scores.txt")
    os.system("sed -i '$d' txt/y_scores.txt")
    
    with open('txt/y_scores.txt', 'r') as file:
        lines = file.readlines()
        for line in lines:
            words = line.split()
            y_angle = words[0].split("_",5)[3]
            y_values.append(float(y_angle))
            scores.append(float(words[1]))
    
    df1 = pd.DataFrame({'scores': scores})
    df2 = pd.DataFrame({'y_angle': y_values})
    df3 = df2.join(df1)
    
    df3 = df3.sort_values(['y_angle'], ascending=True)
    
    plt.scatter(df3['y_angle'], df3['scores'])
    plt.xlabel("Tilt angle (A)")
    plt.ylabel("Scores (REU)")
    plt.savefig("y-angle_plot.png", dpi=600)
    #plt.show()
    
    best_df3 = min(df3['scores'])
    idx = df3.loc[df3['scores'] == best_df3].index[0]
    best_y_angle = int(df3.loc[idx][0])
    return best_y_angle 


def AH_CA_rmsd(native, best):
    number_of_residues = native.size()
    pymol = PyMOLMover()
    
    best_clone = best.clone()
    native_clone = native.clone()
    
    #initiate the spanless membrane mover
    fm = AddSpanlessMembraneMover()
    fm.add_membrane_virtual(native_clone)
    fm.apply(native_clone)
    
    #move the protein along the x and y axes to keep embedding the same
    cmass = pyrosetta.rosetta.core.pose.center_of_mass(native_clone, 1, native_clone.size()-1)
    move_xy = pyrosetta.rosetta.numeric.xyzVector_double_t(-cmass[0], -cmass[1], 0)
    
    #run the translation
    shifted = native_clone.clone()
    translation_mover = pyrosetta.rosetta.protocols.rigid.WholeBodyTranslationMover(move_xy)
    translation_mover.apply(shifted)
    
    ht = HelixTools()
    helix_normal = ht.calculate_screw_axis(shifted)
    angle_with_x = ht.calc_angle(helix_normal,'x')
    angle_with_y = ht.calc_angle(helix_normal,'y')
    
    
    if angle_with_y <= 90:
        rot_z  = pyrosetta.rosetta.protocols.rigid.WholeBodyRotationMover(pyrosetta.rosetta.numeric.xyzVector_double_t(0,0,1), 
                                                                          pyrosetta.rosetta.core.pose.center_of_mass(shifted, 1, shifted.size()-1),
                                                                          -angle_with_x)
    else:
        rot_z  = pyrosetta.rosetta.protocols.rigid.WholeBodyRotationMover(pyrosetta.rosetta.numeric.xyzVector_double_t(0,0,1), 
                                                                          pyrosetta.rosetta.core.pose.center_of_mass(shifted, 1, shifted.size()-1),
                                                                          angle_with_x)
    
    rot_z.apply(shifted)
    shifted.dump_pdb("results/{}/output_pdbs/shifted_{}.pdb".format(line, line))
    
    ca_shifted = []
    ca_best = []

    #calculate the CA positions for the shifted original and the best structures
    for i in range(1,shifted.size()):
        ca_shifted.append(np.array(shifted.residue(i).xyz('CA')))
    for j in range(1, best.size()):
        ca_best.append(np.array(best_clone.residue(j).xyz('CA')))

    total = 0
    for k in range(number_of_residues):
        total += np.square((abs(ca_shifted[k][0]) - abs(ca_best[k][0]))) + np.square((abs(ca_shifted[k][1]) - abs(ca_best[k][1]))) + np.square((abs(ca_shifted[k][2]) - abs(ca_best[k][2])))
        
    average = total/len(ca_best)
    squared_average = np.sqrt(average)
    
    return squared_average


def angle_RMSD(value_list, ideal_value):
    """ Simple RMSD calculator from a 1D angle list and an ideal angle value. """
    total = 0
    for k in range(len(value_list)):
        total = total + np.sqrt(np.square(value_list[k] - ideal_value))
    average = total/len(value_list)
    angle_rmsd = np.sqrt(average)
        
    return angle_rmsd


def get_closest_distance_to_membrane(pose):
    """ Calculate the distance between the membrane core and the helix atom
        closest to the membrane core. """
    pose_clone = pose.clone()
    membrane_mover = AddSpanlessMembraneMover()
    membrane_mover.add_membrane_virtual(pose_clone)
    membrane_mover.apply(pose_clone)
    
    closest_distance = 99999
    for resno in range(1, pose_clone.size()):
        for atom in range(1, pose_clone.residue(resno).natoms()+1):
            distance_from_surface = abs(pose_clone.residue(resno).xyz(atom)[2])
            if distance_from_surface < closest_distance:
                closest_distance = distance_from_surface
    
    return closest_distance
    

def get_HM_vector_angle(native_pose):
    """ Calculates the angle between the HM vector of a native helix pose
        and the z-axis of the coordinate system. """
    native_pose_clone = native_pose.clone() 
    membrane_mover = AddSpanlessMembraneMover()
    membrane_mover.apply(native_pose_clone)
    ht = HelixTools()
    hm = HydrophobicMoment()
    hm_vector = hm.calculate_hydrophobic_moment(native_pose_clone)
    hm_vector_angle = ht.calc_angle(hm_vector, 'z')
    
    return hm_vector_angle

rr = []
name = []
def get_helix_interaction_energy(full_pose, native_pose, sfxn=get_fa_scorefxn()):
    """ Takes the native pose and the helix residue range, and calculates
        the interaction energy between the helix and the rest of the protein.
        The residue range numbering must be done based on the original pdb numbering. """    
    #calculate the residue numbers for the overlapping regions
    seq1 = pyrosetta.rosetta.core.sequence.Sequence()
    seq1.sequence(full_pose.sequence())
    seq2 = pyrosetta.rosetta.core.sequence.Sequence()
    seq2.sequence(native_pose.sequence())
    naive_align = pyrosetta.rosetta.core.sequence.align_naive(seq1, seq2)
    first_res = int(str(naive_align.sequence(1)).split()[1])
    last_res = first_res + naive_align.length() - 1
    residue_range = "{}-{}".format(first_res, last_res)
    rr.append(residue_range)
    name.append(native_pose)
    #define the helix and non-helix regions. 
    select_helix = pyrosetta.rosetta.core.select.residue_selector.ResidueIndexSelector()
    select_helix.set_index(residue_range)    
    select_not_helix = pyrosetta.rosetta.core.select.residue_selector.NotResidueSelector()
    select_not_helix.set_residue_selector(select_helix)

    sfxn(full_pose)
    metric = pyrosetta.rosetta.core.simple_metrics.metrics.InteractionEnergyMetric()
    metric.set_residue_selectors(select_helix, select_not_helix)
    interaction_energy = metric.calculate(full_pose)
    
    return interaction_energy
    

def get_residue_range(full_pose, native_pose):
    """ Takes the native pose and the helix residue range, and calculates
        the interaction energy between the helix and the rest of the protein.
        The residue range numbering must be done based on the original pdb numbering. """    
    #calculate the residue numbers for the overlapping regions
    rr = []; name = []
    seq1 = pyrosetta.rosetta.core.sequence.Sequence()
    seq1.sequence(full_pose.sequence())
    seq2 = pyrosetta.rosetta.core.sequence.Sequence()
    seq2.sequence(native_pose.sequence())
    naive_align = pyrosetta.rosetta.core.sequence.align_naive(seq1, seq2)
    first_res = int(str(naive_align.sequence(1)).split()[1])
    last_res = first_res + naive_align.length() - 1
    residue_range = "{}-{}".format(first_res, last_res)
    rr.append(residue_range)
    #makeshift solution to splice the pose info into file name
    name.append(native_pose.pdb_info().name().split("/")[1])
    
    df = pd.DataFrame()
    df['Name'] = name
    df['Residue range'] = rr
    df.to_csv('residue_range.csv', mode='a', header=False, index=False)
    

def get_membrane_score(pose):
    pose_copy = pose.clone()
    #initiate the spanless membrane mover
    fm = AddSpanlessMembraneMover()
    fm.add_membrane_virtual(pose_copy)
    fm.apply(pose_copy)
    
    mem_sfxn = pyrosetta.get_fa_scorefxn()
    mem_sfxn.add_weights_from_file("mpframework_smooth_fa_2012.wts")

    return mem_sfxn(pose_copy)
    
    
###############################################################################
names = []; rmsds = []; charges = []; sizes = []; tilt_angles = []
psi_angle_rmsd = []; phi_angle_rmsd = []; native_distance_to_membrane = []
native_hm_vector_angles = []; best_hm_vector_angles = []; charged_residues = []
positively_charged_residues = []; negatively_charged_residues = []; polar_residues = []
apolar_residues = []; aliphatic_residues = []; aromatic_residues = []

#this needs to go to a better home :3
if os.path.exists('residue_range.csv'):
    os.remove('residue_range.csv')

with open('folder_list.txt', 'r') as flist:
    lines = flist.read().splitlines()
    for line in lines:
        full = "/home/gulseva/AH_project/benchmark/peripheral_testing/input_pdbs/clean/{}_clean.pdb".format(line.split("_",1)[0])
        native = "results/{}/output_pdbs/{}.pdb".format(line, line)
        best = "results/{}/output_pdbs/best_{}_0.pdb".format(line, line)
        native_pose = pose_from_pdb(native)
        best_pose = pose_from_pdb(best)
                           
        #membrane is initiated at this point        
        rmsds.append(AH_CA_rmsd(native_pose, best_pose))                     
        tilt_angles.append(get_tilt_angle(best_pose))

        thickness_from_file = 0
        with open('metric_data/membrane_thickness.txt', 'r') as thickness_file:
            for thicc in thickness_file:
                if line in thicc:
                    thickness_from_file = round(float(thicc.split(',', 1)[1]))

        calculate_embedded_residues(line, thickness_from_file)















