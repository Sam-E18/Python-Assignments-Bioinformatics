######### Assignment_6: PDB Mean Minimum Distance Calculation ############
    ## *This script reads an ATOM record from a PDB file, groups the atomic coordinates*
    ## *by chain and residue, calculates the minimum distance between every unique pair*
    ## *of residues in each chain (the shortest atom-atom distance), and finally*
    ## *returns the average of these minimum distances for each chain.*
    
import math
import sys

def calculate_pdb_chain_mean_minimum_distances(pdb_file_path=None):
    """
    Calculates the mean of the minimum inter-residue distances for each chain in a PDB file.
    Note: Can be called via 'cat file.pdb | python script.py' or 'python script.py file.pdb'
    """
    
    # --- 1. File Handling ---
    # Read from standard input (stdin) if no file path is provided, otherwise open the file
    if pdb_file_path is None:
        input_fd = sys.stdin
    else:
        # It's good practice to wrap file handling in try/except or use 'with open()'
        # For simplicity, we stick to the provided structure.
        input_fd = open(pdb_file_path, "r")

    # pdb_data stores the parsed coordinates: {chain_id: {res_seq_num: [list_of_coord_tuples]}}
    pdb_data = {}  
    
    # --- 2. PDB File Parsing ---
    for line in input_fd:      
        # Only process lines describing atoms in the structure
        if line.startswith("ATOM"):                         
            # Extract Chain ID (column 22)
            chain_id = line[21].strip()                  
            # Extract Residue Sequence Number (columns 23-26)
            # We strip to handle cases with alignment/spaces
            residue_seq = line[22:26].strip()                
            
            # Extract Coordinates (X, Y, Z) - Columns 31-54
            # Store them as a simple tuple of floats for efficiency
            try:
                x_coord = float(line[30:38].strip())
                y_coord = float(line[38:46].strip())
                z_coord = float(line[46:54].strip())
            except ValueError:
                 # Skip lines where coordinate parsing fails (e.g., malformed PDB)
                continue
            
            coordinates = (x_coord, y_coord, z_coord)
            
            # Pythonic way to build nested dictionaries: use .setdefault()
            # If chain_id is new, default value is {} (empty dictionary)
            # If residue_seq is new in that chain, default value is [] (empty list)
            # Then append the coordinate tuple to the list
            pdb_data.setdefault(chain_id, {}).setdefault(residue_seq, []).append(coordinates)

    # Close the file handle if it was opened from a path (not stdin)
    if pdb_file_path is not None:
        input_fd.close()

    # --- 3. Distance Calculation ---
    result_mean_distances = {} # Will store {chain_id: mean_min_distance}
    
    # Iterate over each protein chain
    for chain_id, residues_dict in pdb_data.items():  
        distances_list = []  # List to hold the minimum distance for EVERY unique residue pair
        
        # Get a list of all residue IDs for the current chain
        residues_list = list(residues_dict.keys())   

        # Double loop to iterate over all UNIQUE pairs of residues (R1 and R2)
        # The 'j in range(i)' pattern ensures we calculate Rj-Ri only once
        for i in range(len(residues_list)):
            residue1_id = residues_list[i]
            
            # Loop over indices j that are less than i to ensure we only compare unique pairs
            for j in range(i):
                residue2_id = residues_list[j]
                
                # List to store distances between all atom pairs of the two residues
                atom_pair_distances = []                 
                
                # Get the lists of coordinates for the two residues
                atoms1_coords = residues_dict[residue1_id]
                atoms2_coords = residues_dict[residue2_id]
                
                # Compare every atom of residue 1 with every atom of residue 2
                for x1, y1, z1 in atoms1_coords:
                    for x2, y2, z2 in atoms2_coords:
                        # Euclidean distance calculation
                        distance_sq = (x1 - x2)**2 + (y1 - y2)**2 + (z1 - z2)**2
                        distance = math.sqrt(distance_sq)
                        atom_pair_distances.append(distance)

                # The distance between the two residues is the minimum distance between their atoms
                if atom_pair_distances:
                    distances_list.append(min(atom_pair_distances))

        # --- 4. Final Calculation and Store Result ---
        # Calculate the mean (average) of all the minimum inter-residue distances
        if len(distances_list) > 0:
            mean_distance = sum(distances_list) / len(distances_list)
            # Store the final result in the dictionary
            result_mean_distances[chain_id] = mean_distance

    # The function returns the final dictionary, separating calculation from I/O.
    return result_mean_distances

# --- Main Execution Block ---
# This block ensures the function is only called when the script is run directly (not imported as a module)
if __name__ == "__main__":
    
    # Check if a file path was provided as a command-line argument (sys.argv[1])
    pdb_file_path = sys.argv[1] if len(sys.argv) > 1 else None
    
    # Calculate the mean minimum distances per chain
    results = calculate_pdb_chain_mean_minimum_distances(pdb_file_path)

    # Print the final output in the required format (Chain: X.XXXX)
    # Using the standard C-style string formatting for cleaner output control
    for chain_id, mean_distance in results.items():
        print("%s: %.4f" % (chain_id, mean_distance), file=sys.stdout)


### Results thats you should obtain

# 1. 1a28.pdb
# - A: 20.4524
# - B: 20.3459

# 2. 1vzq.pdb 
# - H: 18.3158
# - I: 5.6829
# - L: 5.6400

# 3. 1ov9.pdb 
# - A: 13.4893
# - B: 14.8209

# 4. 1jds.pdb 
# - A: 18.7709
# - B: 18.7798
# - C: 19.1461
# - D: 18.7932
# - E: 18.7648
# - F: 18.9342

# 5. UBQ.pdb
# - A: 12.3176

