
#######   Assignment_3:    ############

# This script calculates the frequency of specific amino acid subsequences in a given multiline FASTA protein file (short_examp_fastafile.fa). 
# It reads subsequences from a separate file (sequence_fragments.txt), checks how many protein sequences contain each subsequence at least 10 times (number_of_repetitions), and calculates the proportion of proteins meeting this condition. 
# The results are saved in an output file (exercise3_output1.txt), sorted in descending order by proportion. The output file includes:
# The total number of proteins analyzed.
# The total number of subsequences checked.
# A list of subsequences, their occurrence count, and their proportion, formatted with right alignment.

def calculate_aminoacid_frequencies(fasta_filename, subsequences_filename, number_of_repetitions, output_filename):
    """Given a multiline FASTA_prot file (fasta_filename) and a sub-sequences file (subsequences_filename)
    (one sequence_prot in each line), calculates the proportion of proteins in the FASTA file containing at least N times (number_of repetitions)
    each of the subsequences. Output should bea file sorted by proportion in descending order, with a specific format(rightaligned(pos40), the oder aligned at position 20)."""
    
    ##create the dictionary to the subsequences 
    subsequences = {} ##empty dictionary 
    number_of_repetitions = int(number_of_repetitions)

    with open(subsequences_filename , 'r') as subsequences_filename:
        for subseq_prot in subsequences_filename:
            subsequences[subseq_prot.strip()] = 0    ##add subsequences as keys

     ##count subsequences in each sequence_prot
    protein_counter = 0             
    seq_prot = "" 
    with open(fasta_filename , 'r') as fasta_filename:                        
        for line in fasta_filename:                                                
            if line.startswith(">"):
                if seq_prot:                                        
                    protein_counter += 1
                    for subseq_prot in subsequences: 
                        if seq_prot.count(subseq_prot) >= number_of_repetitions:                    
                            subsequences[subseq_prot] += 1
                seq_prot = ""
            else:
                seq_prot += line.strip()  ##join lines into a single sequence string

        ##count the last protein
        if seq_prot:
            protein_counter += 1
            for subseq_prot in subsequences:                
                if seq_prot.count(subseq_prot) >= number_of_repetitions:                    
                    subsequences[subseq_prot] += 1

    ## write the output file (sort subsequences to proportion in descending order
    sorted_subsequences = sorted(
        subsequences.items() ,
        key = lambda item: item[1],
        reverse=True )
    
    ##write the output file with the headers Number of proteins,Number of subsequences,etc.
    with open(output_filename , 'w') as output_file:
        output_file.write("#Number of proteins:" + f"{protein_counter:>21}\n")
        output_file.write("#Number of subsequences:" + f"{len(subsequences):>17}\n")
        output_file.write("#Subsequence proportions:\n")

        for subseq, count in sorted_subsequences:
            proportion = count / protein_counter
            dist = 20 - len(subseq)
            output_file.write(subseq + f"{count:>{dist}}" + f"{proportion:>20.4f}\n")

    print(f"Output written to {output_filename}")
    return output_filename

print(calculate_aminoacid_frequencies("short_examp_fastafile.fa", "sequence_fragments.txt", "10", "exercise3_output1.txt"))


