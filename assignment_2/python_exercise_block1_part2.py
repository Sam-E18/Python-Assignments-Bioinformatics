#######   Assignment_2:    ############
# In the directory there are the files that the  generate with this script ().
# For the first (1) exercise of this assignment, if you run this script, you must give it the .fa file and if you want to change which residue you want to calculate the frequency depending on the treshold assigned to it. 
# The result is a float number that represents the frequency of the residue.

# For the second exercise (Exercise_2: Output1.txt), you must give the .fa file, the name of the outputfile and the positions you want to display (Example: The first 5 and the last 5). 
# In the directory are the files generated with this script and the test file that was used.


#1)#Calculate the ratio of proteins meeting residue tresholds in a FASTA file.#


def get_proteins_ratio_by_residue_threshold(filename, residue, relative_threshold = 0.03, absolute_threshold= 10):
    """Given a multi-line protein FASTA file (stored in a file with path defined filename), returns a
float corresponding to the ratio of proteins in the fasta file having a relative frequency higher
or equal than a given threshold provided as an argument named “relative_threshold” and
having an absolute frequency of the same residue higher or equal than a given threshold
provided as an argument named “absolute_threshold” for a given residue."""

    with open(filename, "r") as input_file:
        prot_sequences = []
        current_p_seq = ""

        for line in input_file:
            if line.startswith(">"):
                if current_p_seq:  # save previous sequence
                    prot_sequences.append(current_p_seq)
                current_p_seq = ""
            else:
                current_p_seq += line.strip()
        if current_p_seq:  # add the last sequence
            prot_sequences.append(current_p_seq)

    match_proteins_count = 0

    for prot in prot_sequences:
        absolute_count = prot.count(residue)
        relative_frequency = absolute_count / len(prot)

        if relative_frequency >= relative_threshold and absolute_count >= absolute_threshold:
            match_proteins_count += 1

    total_proteins = len(prot_sequences)

    return (match_proteins_count / total_proteins)  #return the final ratio

print(get_proteins_ratio_by_residue_threshold("short_examp_fastafile.fa", residue="M"))



#2)############
def print_sequence_summary(filename , output_filename , first_n = 10 , last_m = 10):
    '''Given a protein fasta file , save on a output file named output_filename the protein identifier, 
    the first n aminoacids, the last M-aminoacids and the absolute frequency in the protein of all 
    the aminoacids found in the protein (the aminoacids that do not appear in the protein should not be shown). 
    The fields must be separated by a tabulator, and one protein by line.'''
    
    with open(filename , "r") as input_file:
        with open(output_filename , "w") as p_output_file:
            p_output_file.write("### Protein Summary Output ###\n")
            p_output_file.write("Header\tFirstN \tLastM \tAmino Acid Counts\n")
            p_output_file.write("----------------------------------------------------------------------------------------------------\n")

            p_sequence = ""  # beggin th p_seq
            p_header = ""  #initialize the header of the p_seq

            for line in input_file:     
                
                if line.startswith(">"):     # is a header      
                    if p_sequence:
                        first_residue = p_sequence[:first_n] + "\t"   
                        last_residue = p_sequence[-last_m:] + "\t"    
                        
                        aa_counts = {}
                        for aa in p_sequence:    #count the residues 
                            if aa in aa_counts:
                                aa_counts[aa] += 1
                            else:
                                aa_counts[aa] = 1
                        
                        count = str(aa_counts)
                        count = count[count.index("{")+1:count.index("}")] + "\n"    
                        count = count.replace("\'" ,"")
                        p_output_file.write(p_header+first_residue+last_residue+count)       
                        
                        p_sequence = ""                                                            
                
                    
                    p_header = line.strip()     # show the header of each proten
                    p_header = p_header.replace(">" , "") + "\t" #replace the > to a empty space
                
                else:                                                                      
                    p_sequence += line.strip()                                                     

            first_residue = p_sequence[:first_n] + "\t" 
            last_residue = p_sequence[-last_m:] + "\t"                                          
                        
            aa_counts = {}
            for aa in p_sequence:                     
                if aa in aa_counts:
                    aa_counts[aa] += 1
                            
                else:
                    aa_counts[aa] = 1
                        
            count = str(aa_counts)        
            count = count[count.index("{")+1:count.index("}")] + "\n"      
            count = count.replace("\'" ,"")
            p_output_file.write(p_header+first_residue+last_residue+count) 
            

print_sequence_summary("short_examp_fastafile.fa", "output2.txt", first_n=10, last_m=4)

