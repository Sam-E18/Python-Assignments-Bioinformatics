#######   Assignment_4:    ############

#This Python script analyzes multiple FASTA files to compare protein identifiers across them. 
# It takes a list of FASTA files as input, where each file contains multiple protein sequences with unique identifiers. 
# The script processes these files using a generator function (FASTA_iterator), which extracts the identifiers and sequences. 
# The main function, compare_fasta_file_identifiers, then compares the identifiers across all input files and generates four key outputs:
#  (1) "intersection," which contains identifiers found in all files; (2) "union," which includes all unique identifiers across the files;
#  (3) "frequency," a dictionary that records how many files contain each identifier; 
# and (4) "specific," a dictionary mapping each FASTA file to a set of identifiers that are unique to that file.
#  The results are written to an output file named pytexercise_4_output.txt, where they are structured for easy interpretation. 
# The expected output includes sets of identifiers and a frequency count, allowing users to understand the distribution of protein identifiers across the FASTA files.


## 1) A Generator Function that reads a Fasta file. 
## Each iteration returns a tuple with the format: (identifier, sequence).

def FASTA_iterator(fasta_filename):
    """Generator function that reads a FASTA file line by line.
    
    Each iteration returns a tuple (identifier, sequence).
    The identifier is extracted from lines starting with '>',
    and the sequence is concatenated from subsequent lines.
    """
    
    current_id = ""  # Stores the identifier of a sequence
    current_p_seq = ""  # Stores the sequence itself

    with open(fasta_filename, "r") as input_file:  # Open file safely
        for line in input_file:
            if line.startswith(">"):  # If a new identifier is found
                if current_p_seq:  
                    yield (current_id, current_p_seq)  # Yield previous sequence before resetting
                    current_p_seq = ""  # Reset sequence storage
                current_id = line.strip().replace(">", "")  # Extract the identifier (without '>')
            else:       
                current_p_seq += line.strip()  # Append sequence lines

        yield (current_id, current_p_seq)  # Yield the last sequence in the file

## 2) Compare FASTA file identifiers and return a dictionary containing 4 key results.

def compare_fasta_file_identifiers(fasta_filenames_list):
    """Given a list of FASTA files, this function returns a dictionary containing:
    
    - "intersection": A set of common identifiers found in all FASTA files.
    - "union": A set of all unique identifiers found across all FASTA files.
    - "frequency": A dictionary with identifiers as keys and the number of files they appear in as values.
    - "specific": A dictionary with filenames as keys and sets of identifiers unique to that file.
    """

    ## Dictionaries to store results
    file_to_ids = {}  # Maps each file to a set of its identifiers
    frequency_dict_ids = {}  # Tracks the frequency of each identifier across files
    specific_dict_ids = {}  # Stores unique identifiers per file
    output_dictionary = {}  # Final dictionary to return results

    # Step 1: Read all FASTA files and extract identifiers
    for file in fasta_filenames_list:     
        set_ids_file = set()  # Stores unique identifiers for this file
        
        for id_protein, _ in FASTA_iterator(file):  # Iterate through sequences in the file
            set_ids_file.add(id_protein.upper())  # Convert to uppercase for case insensitivity
            if id_protein.upper() not in frequency_dict_ids:
                frequency_dict_ids[id_protein.upper()] = 1  # First occurrence of identifier
            else:
                frequency_dict_ids[id_protein.upper()] += 1  # Increase frequency count

        file_to_ids[file] = set_ids_file  # Store extracted IDs in the dictionary
        specific_dict_ids[file] = set()  # Initialize dictionary for specific identifiers

    list_of_id_sets = list(file_to_ids.values())  # Convert dictionary values to a list of sets

    # Step 2: Detect file-specific identifiers
    for index, file_name in enumerate(file_to_ids):
        # Create a union of identifiers from all other files except the current one
        other_files_union = set().union(*[s for i, s in enumerate(list_of_id_sets) if i != index])
        specific_dict_ids[file_name] = file_to_ids[file_name].difference(other_files_union)  # Find unique IDs

    # Step 3: Store results in the output dictionary
    output_dictionary["intersection"] = set.intersection(*list_of_id_sets)  # Common identifiers in all files
    output_dictionary["union"] = set.union(*list_of_id_sets)  # All unique identifiers
    output_dictionary["frequency"] = frequency_dict_ids  # How many files contain each identifier
    output_dictionary["specific"] = specific_dict_ids  # Unique identifiers per file
    
    return output_dictionary  # Return the result dictionary

# List of FASTA files to analyze
fasta_files = [
    "1_uniprot_sprot_short.fasta",
    "2_uniprot_sprot_short.fasta",
    "3_uniprot_sprot_short.fasta"
]

# Compute identifier comparison
result = compare_fasta_file_identifiers(fasta_files)

# Write results to a file
output_filename = "pytexercise_4_output.txt"

with open(output_filename, "w") as output_file:
    output_file.write("🔹 INTERSECTION (Common Identifiers in All Files):\n")
    output_file.write(f"{result['intersection']}\n\n")

    output_file.write("🔹 UNION (All Unique Identifiers Across Files):\n")
    output_file.write(f"{result['union']}\n\n")

    output_file.write("🔹 FREQUENCY (Number of Files Each Identifier Appears In):\n")
    for key, value in result["frequency"].items():
        output_file.write(f"{key}: {value}\n")

    output_file.write("\n🔹 SPECIFIC (Unique Identifiers Per File):\n")
    for filename, unique_ids in result["specific"].items():
        output_file.write(f"{filename}: {unique_ids}\n")

# Print confirmation message
print(f"Results have been written to {output_filename}")
