######### Assignment_6: PDB Mean Minimum Distance Calculation ############
    ## *This script reads an ATOM record from a PDB file, groups the atomic coordinates*
    ## *by chain and residue, calculates the minimum distance between every unique pair*
    ## *of residues in each chain (the shortest atom-atom distance), and finally*
    ## *returns the average of these minimum distances for each chain.*
    
import math
import sys

def calculate_pdb_chain_mean_minimum_distances(pdbfile_path=None):
    """Calculate the mean of the minimum distance between any two residues pairs found in the same chain of a PDB.
    Note: call this script from the terminal with:1. cat file.pdb | python script.py  2. python script.py file.pdb"""
       # call the file from the stdin
    if pdbfile_path is None:
        structure_pdb = sys.stdin
    else:
        structure_pdb = open(pdbfile_path, "r")

    output_dictionary = {} 
    # take the specific data from the colums and put in dictionaries
    for line in structure_pdb:      
        if line.startswith("ATOM"):                         
            atom_name = line[13:16].strip()               
            residue_prot = line[17:20].strip()               
            chain_id = line[21].strip()                  
            residue_seq = line[22:26].strip()                
            coords = {'X': float(line[30:38].strip()),   
                     'Y': float(line[38:46].strip()),
                     'Z': float(line[46:54].strip())}
            
            if chain_id not in output_dictionary:
                chain_dictionary = {}
                output_dictionary[chain_id] = chain_dictionary

            if residue_prot + residue_seq not in chain_dictionary:
                residue_dictionary = {}

            residue_dictionary[atom_name] = coords            
            chain_dictionary[residue_prot+residue_seq] = residue_dictionary  
            output_dictionary[chain_id] = chain_dictionary
       #iterate with for loops around the residues
    chain_distance = {} 
    for chain_id in output_dictionary:  
        residues_distance = {}  
        residues_list = list(output_dictionary[chain_id].keys())   

        for i in range(len(residues_list)):
            residue1 = residues_list[i]

            for residue2 in residues_list[i+1:]:
                atoms_distance = []                 
                for atom1 in output_dictionary[chain_id][residue1]:
                    for atom2 in output_dictionary[chain_id][residue2]:
                        coordinates1 = output_dictionary[chain_id][residue1][atom1]
                        coordinates2 = output_dictionary[chain_id][residue2][atom2]

                        atoms_distance.append(math.sqrt((coordinates1['X']-coordinates2['X'])**2 +
                                                       (coordinates1['Y']-coordinates2['Y'])**2 +
                                                       (coordinates1['Z']-coordinates2['Z'])**2))

                residues_distance[residue1 + "_" + residue2] = min(atoms_distance)

        chain_distance[chain_id] = residues_distance

    for chain_id in chain_distance:
        sum_dist = 0
        count_pairs = 0

        for pair in chain_distance[chain_id]:
            count_pairs += 1
            sum_dist += chain_distance[chain_id][pair]

        chain_mean_distance = sum_dist / count_pairs
        chain_distance[chain_id] = float(f"{chain_mean_distance:.4f}")

    for item in chain_distance.items():
        output_string = str(item)
        output_string = output_string.replace("(", "", -1)
        output_string = output_string.replace(")", "", -1)
        output_string = output_string.replace("'", "", -1)
        output_string = output_string.replace(",", ":", -1)
        print(output_string, file=sys.stdout)

   
    if pdbfile_path is not None:
        structure_pdb.close()

    return 

# run the script independtly and process the pdbfile
if __name__ == "__main__":

    pdb_file_path = sys.argv[1] if len(sys.argv) > 1 else None
    calculate_pdb_chain_mean_minimum_distances(pdb_file_path)
