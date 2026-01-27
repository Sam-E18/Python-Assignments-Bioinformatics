import sequence_dictionaries 
import sys
import os


class IncorrectSequenceLetter(ValueError):
    """New ValueError exception subclass:o create a new exception instance, it should be created with the letter not found in the 
        alphabet and the class name of the sequence.e = IncorrectSequenceLetter(“B”, class_name) """
    def __init__(self, letter, class_type):
        self.letter = letter
        self.class_type = class_type.__name__  ## define the class.

    def __str__(self):
        return f"The sequence item {self.letter} is not found in the alphabet of class {self.class_type}"


def FASTA_iterator(fasta_filename, sequence_class ):
    """ Modifies the FASTA_iterator to generate instances of any subclass of Sequence. 
    Skips sequences with incorrect letters by handling IncorrectSequenceLetter exceptions. """
    with open(fasta_filename, "r") as file_fa:
        any_sequence = ""
        sequence_header = None  
        for line in file_fa:
            if line.startswith(">"): 
                if any_sequence:  
                    try:
                        yield sequence_class(seq_header, any_sequence)
                    except IncorrectSequenceLetter as e:
                        sys.stderr.write(f"{e}\n")
                    any_sequence = ""
                seq_header = line.strip().replace(">", "")
            else:
                any_sequence += line.strip()
        if any_sequence:
            try:
                yield sequence_class(seq_header, any_sequence)
            except IncorrectSequenceLetter as e:
                sys.stderr.write(f"{e}\n")


class Sequence():
    """ Define the differents behaviour for the classes defined ."""
    alphabet = ""
    weights = {}
    def __init__(self, identifier, sequence ):
        self.__identifier = identifier
        self.__sequence = sequence

        for residue_prot_dna_rna in self.__sequence: ## any of the residues form prot dna and rna
            if residue_prot_dna_rna not in self.alphabet:
                raise IncorrectSequenceLetter(residue_prot_dna_rna, self.__class__) ## show the error of incorrect sequence letter.
                
    def get_identifier(self):
        return self.__identifier
    
    def get_sequence(self):
        return self.__sequence 

    
    def get_mw(self): ## its better to put the calculation in one line with the generator.
        return round(sum(self.weights.get(residue, 0) for residue in self.__sequence) - (len(self.__sequence) - 1) * 18.015, 2)

    
    def has_subsequence(self, sequence_query_instance): ## this could be class Sequence or ProteinSequence or DNASequence or RNAsequence.
        return sequence_query_instance.get_sequence() in self.get_sequence()
    
    ## should return the lenght of the sequence len(sequence).
    def __len__(self):
        return len(self.get_sequence())
    
    ## return True if sequence strings are exactly the same (without taking into account the identifiers). 
    ## sequence1 == sequence2: 
    def __eq__(self,other):
        return self.get_sequence() == other.get_sequence()
    
    ## return True if sequences are different, without taking into account the identifiers.
    ## sequence1 != sequence2:
    def __ne__(self, other):
        return self.get_sequence() != other.get_sequence()
    
    ## return a concatenate sequence 
    ## Sequence  +  Sequence: 
    def __add__(self, other):
        if type(self) != type(other):
            raise TypeError("Sequences must be instances of the same class to be added.")
        
        concatenate_identifier = self.get_identifier() + "+" + other.get_identifier()
        concatenate_sequence = self.get_sequence() + other.get_sequence()
        return self.__class__(concatenate_identifier, concatenate_sequence)

    # its necesary to print the concatenate sequence 
    def __str__(self) -> str:
        return self.get_identifier() + " == " + self.get_sequence() ## “\n” is printed below, with the == written next to it.
    
    ##Should  return  the  sequence  element  at  position  i.  Position  0 corresponds to the first position
    ## return a sequence i 
    def __getitem__(self,key):
        return self.get_sequence()[key]
    
    ##in operator: should return a boolean if the string is a substring of the attribute sequence.
    def __contains__(self , item):
        return item.get_sequence() in self.get_sequence()
    
    ##comparing sequences
    """ Implement  the  necessary  method(s)  to  define  how 
    sequences  should  be  ordered.  The  objective  is  that  when  sorting  a  list  of 
    sequences, they are sorted according to their molecular weight."""

    def __lt__(self, other):
        return self.get_mw() < other.get_mw()
    
    def __le__(self,other):

        return self.get_mw() <= other.get_mw()
    
    def __gt__(self,other):

        return self.get_mw() > other.get_mw()
    
    def __ge__(self,other):

        return self.get_mw() >= other.get_mw()
    
    """Adapt the sequence class so that it can be used as key in a dictionary or it can be 
    added to a set. Two sequences should be considered the same object in terms of 
    set or key if they share both the identifier and the sequence."""
    ## Using __hash__
    def __eq__(self, other):
        return self.get_identifier() == other.get_identifier() and self.get_sequence() == other.get_sequence()

    def __hash__(self):
        return hash((self.get_identifier(), self.get_sequence()))
    


class ProteinSequence(Sequence):
    alphabet = sequence_dictionaries.protein_letters
    weights = sequence_dictionaries.protein_weights
   
    

class NucleotideSequence(Sequence):

    translate_table = {}
    stop_codons = []

    def translate(self):
        translated_sequence = ""
        for codon_start in range(0, len(self.get_sequence()), 3):
            codon = self.get_sequence()[codon_start:codon_start + 3]
            if len(codon) < 3:
                continue
            if codon in self.stop_codons:
                translated_sequence += " "  # Standard for stop codons
            elif codon in self.translate_table:
                translated_sequence += self.translate_table[codon]
        return ProteinSequence(self.get_identifier(), translated_sequence or "X")
    
class DNASequence(NucleotideSequence):
    translate_table = sequence_dictionaries.dna_table
    alphabet = sequence_dictionaries.dna_letters
    weights = sequence_dictionaries.dna_weights
    stop_codons = sequence_dictionaries.dna_stop_codons

    def transcribe(self):
        transcribed_seq = self.get_sequence().replace("T", "U")
        return RNASequence(self.get_identifier(), transcribed_seq)


class RNASequence(NucleotideSequence):
    
    translate_table = sequence_dictionaries.rna_table 
    alphabet = sequence_dictionaries.rna_letters 
    weights = sequence_dictionaries.rna_weights
    stop_codons = sequence_dictionaries.rna_stop_codons
    



## MAIN block
""" 7) When the script is executed as a standalone application (without being imported (i.e. 
code under __main__ block), the script should read input DNA FASTA file(s) and calculate 
the length and molecular weight of their corresponding proteins (i.e. corresponding to 
ProteinSequence instances obtained after translation). The script must print the output 
to standard output or to a file. Output should be sorted by molecular weight, from 
lowest to greatest. """

def list_fasta_files(current_directory):
    """7.1) Lists all .fasta and .fa files in the specified directory."""
    fasta_files = []
    for filename in os.listdir(current_directory):
        if filename.endswith(".fasta") or filename.endswith(".fa"):
            fasta_files.append(os.path.join(current_directory, filename))
    return fasta_files

def process_files(input_path, output_file=None):
    """7.2 and 7.3 Processes the input path and writes sorted results to standard output or a file."""
    if os.path.isdir(input_path):
        fasta_files = list_fasta_files(input_path)
    else:
        fasta_files = [input_path] if input_path.endswith((".fasta", ".fa")) else []

    sequences = []
    sys.stderr.write(f"{len(fasta_files)} FASTA files found.\n")

    for fasta_file in fasta_files:
        sys.stderr.write(f"Processing {fasta_file}...\n")
        for sequence in FASTA_iterator(fasta_file, DNASequence):
            try:
                if isinstance(sequence, (DNASequence, RNASequence)):
                    sequence = sequence.translate()  # Ensure it is a ProteinSequence
                sequences.append(sequence)
            except IncorrectSequenceLetter as e:
                sys.stderr.write(f"Skipping sequence {sequence.get_identifier()}: {e}\n")
        sys.stderr.write(f"{fasta_file} finished.\n")
    
    sys.stderr.write(f"{len(sequences)} sequences found.\n")
    sys.stderr.write("Sorting the sequences...Thanks for your patience.\n")
    sequences.sort()
    sys.stderr.write("Sort process finished.*--*\n")

    output = "\n".join(f"{seq.get_identifier()}\t{len(seq)}\t{seq.get_mw()}" for seq in sequences)
    
    if output_file:
        with open(output_file, "w") as f:
            f.write(output + "\n")
    else:
        print(output)
    
    sys.stderr.write("Program finished correctly.(+___+). \n")

if __name__ == "__main__":
    if len(sys.argv) == 1:
        process_files(".")
    elif len(sys.argv) == 2:
        process_files(sys.argv[1])
    elif len(sys.argv) == 3:
        process_files(sys.argv[1], sys.argv[2])
    else:
        sys.stderr.write("Use like: python python_exercise_block2_part4.py [IN] [OUT]\n")
##if you do not supply the input file and outputfile, the script will look for .fasta or .fa files in the current directory and run the program and print it to the terminal.


