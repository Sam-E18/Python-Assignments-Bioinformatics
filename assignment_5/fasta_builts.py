
#uses in the others script: from fasta_builts import FASTA_iterator
def FASTA_iterator(fasta_filename):
    '''Generator function that reads a fasta file. In each iteration, the function
    returns a tuple with the following format: (identifier , sequence).'''
    current_id = ""
    current_p_seq = ""
    with open(fasta_filename, "r") as input_file:
        for line in input_file:
            if line.startswith(">"):                           
                if current_p_seq:
                    yield (current_id, current_p_seq)  # generator with yield
                    current_p_seq = ""                                   
                current_id = line.strip().replace(">", "")
            else:       
                current_p_seq += line.strip()                           
        yield (current_id, current_p_seq)  # yield the last sequence


