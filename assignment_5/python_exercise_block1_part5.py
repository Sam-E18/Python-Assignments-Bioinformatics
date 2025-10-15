
#######   Assignment_5:    ############
## 8 exercise related to each one.

##1)##  Repeat the same exercises proposed in session 2 but using the FASTA_Iterator function created in session 4 to read the FASTA files.

from fasta_builts import FASTA_iterator


# def FASTA_iterator(fasta_filename):
#     """Generator function that reads a fasta file. In each iteration, the function
#     returns a tuple with the following format: (identifier , sequence)."""
#     current_id = ""
#     current_p_seq = ""
#     with open(fasta_filename, "r") as input_file:
#         for line in input_file:
#             if line.startswith(">"):                           
#                 if current_p_seq:
#                     yield (current_id, current_p_seq)  # generator with yield
#                     current_p_seq = ""                                   
#                 current_id = line.strip().replace(">", "")
#             else:       
#                 current_p_seq += line.strip()                           
#         yield (current_id, current_p_seq)  # yield the last sequence

### EXERCISES of the SESSION 2

def get_proteins_ratio_by_residue_threshold(fasta_filename, residue, relative_threshold=0.03, absolute_threshold=10):
    """Given a multi-line protein FASTA file (stored in a file with path defined filename), returns a
    float corresponding to the ratio of proteins in the fasta file having a relative frequency higher
    or equal than a given threshold provided as an argument named “relative_threshold” and
    having an absolute frequency of the same residue higher or equal than a given threshold
    provided as an argument named “absolute_threshold” for a given residue."""
    
    match_proteins_count = 0
    total_proteins = 0

    for identifier, sequence in FASTA_iterator(fasta_filename):
        total_proteins += 1
        absolute_count = sequence.count(residue)
        relative_frequency = absolute_count / len(sequence)

        if relative_frequency >= relative_threshold and absolute_count >= absolute_threshold:
            match_proteins_count += 1

    return match_proteins_count / total_proteins if total_proteins > 0 else 0.0


print(get_proteins_ratio_by_residue_threshold("short_examp_fastafile.fa", residue="M")) ### I try too with this other fileexample_fasta_file_alg.fa, to make sure that it work well.



def print_sequence_summary(fasta_filename, output_filename, first_n=10, last_m=10):
    """Given a protein fasta file, save on an output file named output_filename the protein identifier, 
    the first n amino acids, the last m amino acids, and the absolute frequency in the protein of all 
    the amino acids found in the protein (the amino acids that do not appear in the protein should not be shown). 
    The fields must be separated by a tabulator, and one protein by line."""
    
    with open(output_filename, "w") as p_output_file:
        p_output_file.write("### Protein Summary Output ###\n")
        p_output_file.write("Header\tFirstN \tLastM \tAmino Acid Counts\n")
        p_output_file.write("----------------------------------------------------------------------------------------------------\n")
        
        # take the FASTA_iterator to loop through the protein data
        for p_header, p_sequence in FASTA_iterator(fasta_filename):
            
            first_residue = p_sequence[:first_n] + "\t"
            last_residue = p_sequence[-last_m:] + "\t"
            
            
            aa_counts = {}
            for aa in p_sequence:
                if aa in aa_counts:
                    aa_counts[aa] += 1
                else:
                    aa_counts[aa] = 1

            count = str(aa_counts)
            count = count[count.index("{") + 1:count.index("}")] + "\n"    
            count = count.replace("'", "")  
            
            p_output_file.write(p_header + "\t" + first_residue + last_residue + count)

print_sequence_summary("short_examp_fastafile.fa", "exercise_5_summary.txt", first_n=5, last_m=5)



##2)## A function that, given a multiline FASTA file, returns the length of the sequence with the maximum length

def get_max_sequence_lenght_from_FASTA_file(fasta_filename):
    """A function that, given a multiline FASTA ﬁle, returns the length of the sequence with the maximum length"""
    lenght_dictionary = {}
    for seq_prot in FASTA_iterator(fasta_filename):
        lenght_dictionary[seq_prot[0]] = len(seq_prot[1])

    return(max(lenght_dictionary.values()))

get_max_sequence_lenght_from_FASTA_file("short_examp_fastafile.fa")



##3)## A function that, given a multiline FASTA file, returns the length of the sequence with the minimum length

def get_min_sequence_lenght_from_FASTA_file(fasta_filename):
    """ A function that, given a multiline FASTA ﬁle, returns the length of the sequence with the minimum length"""
    lenght_dictionary = {}
    for seq_prot in FASTA_iterator(fasta_filename):
        lenght_dictionary[seq_prot[0]] = len(seq_prot[1])

    return(min(lenght_dictionary.values()))

get_min_sequence_lenght_from_FASTA_file("short_examp_fastafile.fa")

# ##4)## A function that, given a FASTA file, returns a list of tuples (identifier, sequence) 
# corresponding to the sequence(s) with maximum length. The list must be sorted 
# by the identifier (case insensitive sorted)

def get_longest_sequences_from_FASTA_file( fasta_filename ):
    """A function that, given a FASTA ﬁle, returns a list of tuples (identiﬁer, sequence) 
    corresponding to the sequence(s) with maximum length. The list must be sorted 
    by the identiﬁer (case insensitive sorted)."""
    maximum_length = get_max_sequence_lenght_from_FASTA_file(fasta_filename)

    my_tupple_list = [seq_prot for seq_prot in FASTA_iterator(fasta_filename) if len(seq_prot[1]) == maximum_length]
    my_tupple_list.sort(key = lambda seq_prot: seq_prot[0].lower())
    
    return(my_tupple_list)

get_longest_sequences_from_FASTA_file("short_examp_fastafile.fa")

# ##5)##  function that, given a FASTA file, returns a list of tuples (identifier, sequence) 
# corresponding to the sequence(s) with minimum length. The list must be sorted by 
# the identifier (case insensitive sorted)

def get_shortest_sequences_from_FASTA_file( fasta_filename ):
    """A function that, given a FASTA ﬁle, returns a list of tuples (identiﬁer, sequence) 
    corresponding to the sequence(s) with minimum length. The list must be sorted by 
    the identiﬁer (case insensitive sorted)."""
    minimum_length = get_min_sequence_lenght_from_FASTA_file(fasta_filename)
    my_tupple_list = [seq_prot for seq_prot in FASTA_iterator(fasta_filename) if len(seq_prot[1]) == minimum_length]
    my_tupple_list.sort(key = lambda seq_prot: seq_prot[0].lower())
    return(my_tupple_list)

get_shortest_sequences_from_FASTA_file("short_examp_fastafile.fa")

# ##6)##A function that, given a protein FASTA file, returns a dictionary with the molecular 
# weights  of  all  the  proteins  in  the  file.  The  dictionary  keys  must  be  the  protein 
# identifiers and the associated values must be a float corresponding to the molecular 
# weight.

def get_molecular_weights(fasta_filename):
    """A function that, given a protein FASTA ﬁle, returns a dictionary with the molecular 
    weights of all the proteins in the ﬁle. The dictionary keys must be the protein 
    identiﬁers and the associated values must be a ﬂoat corresponding to the molecular 
    weight."""
    aa_weights = {'A': 89.09,'R': 174.20,'N': 132.12,'D': 133.10,'C': 121.16,'Q': 146.15,'E': 147.13,
              'G': 75.07,'H': 155.16,'I': 131.18,'L': 131.18,'K': 146.19,'M': 149.21,'F': 165.19,
              'P': 115.13,'S': 105.09,'T': 119.12,'W': 204.23,'Y': 181.19,'V': 117.15 }

    protein_dictionary = {}
    for seq_protein in FASTA_iterator(fasta_filename):
        protein_dictionary[seq_protein[0]] = 0 
        for resid in seq_protein[1]:
            if aa_weights.get(resid):
                protein_dictionary[seq_protein[0]] += aa_weights.get(resid)
            else:
                print(f"Error: The resid {resid} is not a aminoacid :(")

    return(protein_dictionary)

print(get_molecular_weights("short_examp_fastafile.fa"))



# ##7)##A  function  that,  given  a  protein  FASTA  file,  returns  a  tuple  with  (identifier, 
# sequence) of the protein with the lowest molecular weight. If there are two or more 
# proteins having the minimum molecular weight, just return the first one.

def get_sequence_with_min_molecular_weight( fasta_filename ):
    """ A function that, given a protein FASTA ﬁle, returns a tuple with (identiﬁer, 
    sequence) of the protein with the lowest molecular weight. If there are two or more 
    proteins having the minimum molecular weight, just return the ﬁrst one."""
    protein_dictionary_weights = get_molecular_weights(fasta_filename)
    minimum_weight = min(protein_dictionary_weights.values())
    for item in protein_dictionary_weights.items():
        if abs(item[1] - minimum_weight)  < 0.001:
            for seq_protein in FASTA_iterator(fasta_filename):
                if seq_protein[0] == item[0]:
                    return_tuple = seq_protein
    return(return_tuple)

print(get_sequence_with_min_molecular_weight("short_examp_fastafile.fa"))

# ##8)##A function that, given a protein FASTA file, returns the mean of the molecular 
# weights of all the proteins.

def get_mean_molecular_weight( fasta_filename ):
    """ A function that, given a protein FASTA ﬁle, returns the mean of the molecular weights of all the proteins"""
    
    protein_dictionary_weights = get_molecular_weights(fasta_filename)
    sum_of_weights = 0 
    for weight in protein_dictionary_weights.values():
        sum_of_weights += weight
    
    mean_weight = sum_of_weights/len(protein_dictionary_weights)
    return (mean_weight)

get_mean_molecular_weight("short_examp_fastafile.fa")


##PD: do a merge  with the other in the notes_5_jupyter notebook