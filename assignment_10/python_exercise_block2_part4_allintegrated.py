######### Assignment_10: Integrated Bioinformatics Pipeline ############
## * This script implements a complete bioinformatics tool for sequence *
## * analysis. It combines OOP, custom exceptions, and biological logic *
## * to process DNA/RNA/Protein FASTA files, performing transcription, *
## * translation (with Start/Stop detection), and MW calculations. *
## * The script handles errors gracefully and provides a user-friendly  *
## * command-line interface using argparse.                          *
##### Some parts of the code are adapted from previous assignments. #####
##### Big steps #####
## 1) Efficiency: In Python, checking if a letter is in a set is much faster than checking a list or string. This is crucial when processing large genome files.
## 2) Biological Logic (Start/Stop Codons): The original translation logic translated the entire sequence. The newlogic is more biologically accurate: it waits for a Start Codon (like ATG) to begin and stops at a Stop Codon. 
## 3) Command Line Interface (argparse): Instead of using sys.argv (which is hard to manage), we now use argparse. This allows us to use professional "flags" like -i for input and -o for output.
## 4) Caching (self.__mw): We now "save" the result of the molecular weight calculation. If we ask for the weight twice, the script doesn't recalculate it, it just gives you the saved number.
##############################################################################


import sequence_dictionaries 
import sys
import os
import argparse
import re

class IncorrectSequenceLetter(ValueError):
    """New ValueError exception subclass to handle invalid biological residues."""
    def __init__(self, letter, class_type):
        self.letter = letter
        # We handle both string and class type input for flexibility
        self.class_name = class_type if isinstance(class_type, str) else class_type.__name__

    def __str__(self):
        return f"The sequence item {self.letter} is not found in the alphabet of class {self.class_name}"


def FASTA_iterator(fasta_filename, sequence_class):
    """
    Modifies the FASTA_iterator to generate instances of any subclass of Sequence. 
    Skips sequences with incorrect letters by handling IncorrectSequenceLetter exceptions.
    """
    with open(fasta_filename, "r") as file_fa:
        any_sequence = ""
        seq_header = None  
        for line in file_fa:
            if line.startswith(">"): 
                if any_sequence:  
                    try:
                        yield sequence_class(seq_header, any_sequence)
                    except IncorrectSequenceLetter as e:
                        sys.stderr.write(f"Error in {seq_header}: {e}\n")
                    any_sequence = ""
                seq_header = line.strip().replace(">", "")
            else:
                any_sequence += line.strip()
        if any_sequence:
            try:
                yield sequence_class(seq_header, any_sequence)
            except IncorrectSequenceLetter as e:
                sys.stderr.write(f"Error in {seq_header}: {e}\n")


class Sequence():
    """Define the behavior for all sequence types (DNA, RNA, Protein)."""
    alphabet = set() # Use set for O(1) lookup speed
    weights = {}
    
    def __init__(self, identifier, sequence):
        self.__identifier = identifier
        self.__sequence = sequence.upper()
        self.__mw = None # Cache for molecular weight

        # Validation logic using professor's set approach
        for residue_prot_dna_rna in self.__sequence:
            if residue_prot_dna_rna not in self.alphabet:
                raise IncorrectSequenceLetter(residue_prot_dna_rna, self.__class__)
                
    def get_identifier(self):
        return self.__identifier
    
    def get_sequence(self):
        return self.__sequence 

    def get_mw(self):
        """Calculates molecular weight. Base class handles nucleotides."""
        if self.__mw is None:
            self.__mw = sum(self.weights.get(res, 0) for res in self.__sequence)
        return self.__mw

    def __len__(self):
        return len(self.get_sequence())
    
    def __str__(self) -> str:
        return f"{self.get_identifier()} == {self.get_sequence()}"

    def __add__(self, other):
        if type(self) != type(other):
            raise TypeError("Only sequences of the same type can be concatenated.")
        return self.__class__(f"{self.get_identifier()}+{other.get_identifier()}", 
                              self.get_sequence() + other.get_sequence())

    def __lt__(self, other):
        """Sorts sequences by Molecular Weight (as per your requirement)."""
        return self.get_mw() < other.get_mw()

    def __eq__(self, other):
        return self.get_identifier() == other.get_identifier() and self.get_sequence() == other.get_sequence()

    def __hash__(self):
        return hash((self.get_identifier(), self.get_sequence()))


class ProteinSequence(Sequence):
    """Protein subclass with water loss correction for peptide bonds."""
    alphabet = set(sequence_dictionaries.protein_letters)
    weights = sequence_dictionaries.protein_weights

    def get_mw(self):
        base_mw = super().get_mw()
        if len(self) > 1:
            return round(base_mw - (len(self) - 1) * 18.015, 2)
        return round(base_mw, 2)


class NucleotideSequence(Sequence):
    """Base class for DNA and RNA translation logic."""
    translate_table = {}
    stop_codons = set()
    start_codons = set()

    def translate(self):
        """Biologically accurate translation: Start Codon -> Stop Codon."""
        seq = self.get_sequence()
        translated_str = ""
        started = False
        
        for i in range(0, len(seq) - (len(seq) % 3), 3):
            codon = seq[i:i+3]
            if not started:
                if codon in self.start_codons:
                    started = True
                    translated_str += self.translate_table.get(codon, "")
            else:
                if codon in self.stop_codons:
                    break
                translated_str += self.translate_table.get(codon, "")
        
        return ProteinSequence(self.get_identifier(), translated_str or "X")


class DNASequence(NucleotideSequence):
    alphabet = set(sequence_dictionaries.dna_letters)
    weights = sequence_dictionaries.dna_weights
    translate_table = sequence_dictionaries.dna_table
    stop_codons = set(sequence_dictionaries.dna_stop_codons)
    start_codons = set(sequence_dictionaries.dna_start_codons)

    def transcribe(self):
        return RNASequence(self.get_identifier(), self.get_sequence().replace("T", "U"))


class RNASequence(NucleotideSequence):
    alphabet = set(sequence_dictionaries.rna_letters)
    weights = sequence_dictionaries.rna_weights
    translate_table = sequence_dictionaries.rna_table
    stop_codons = set(sequence_dictionaries.rna_stop_codons)
    start_codons = set(sequence_dictionaries.rna_start_codons)


def main():
    # Setting up ARGPARSE for command-line interface
    parser = argparse.ArgumentParser(description="Bioinformatics tool to process FASTA files and calculate Protein MW.")
    parser.add_argument('-i', '--input', dest="infile", default="./", help="Input FASTA file or directory")
    parser.add_argument('-o', '--output', dest="outfile", default=None, help="Output file (optional)")
    parser.add_argument('-v', '--verbose', action="store_true", help="Print progression to stderr")
    
    args = parser.parse_args()

    # Identifying files
    if os.path.isdir(args.infile):
        files = [os.path.join(args.infile, f) for f in os.listdir(args.infile) if f.endswith((".fa", ".fasta"))]
    else:
        files = [args.infile]

    if args.verbose:
        sys.stderr.write(f"Found {len(files)} files.\n")

    results = []
    for f in files:
        for dna in FASTA_iterator(f, DNASequence):
            # We convert DNA -> Protein directly
            protein = dna.translate()
            results.append(protein)
        if args.verbose: sys.stderr.write(f"Finished processing {f}\n")

    # Sorting by MW
    results.sort()

    # Output handling
    output_data = "\n".join(f"{p.get_identifier()}\t{len(p)}\t{p.get_mw()} g/mol" for p in results)
    
    if args.outfile:
        with open(args.outfile, "w") as out:
            out.write(output_data + "\n")
    else:
        print(output_data)

if __name__ == "__main__":
    main()