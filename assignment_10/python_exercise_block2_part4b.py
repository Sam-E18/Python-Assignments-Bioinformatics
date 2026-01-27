import sequence_dictionaries
import sys
import os

class IncorrectSequenceLetter(ValueError):
    """Custom exception for sequences with incorrect letters."""
    def __init__(self, letter, class_type):
        self.letter = letter
        self.class_type = class_type.__name__  # Corrected to get the actual class name

    def __str__(self):
        return f"The sequence item {self.letter} is not found in the alphabet of class {self.class_type}"

def FASTA_iterator(fasta_filename):
    """Detects and yields the correct type of sequence: DNASequence, RNASequence, or ProteinSequence."""
    with open(fasta_filename, "r") as file_fa:
        any_sequence = ""
        seq_header = None
        for line in file_fa:
            if line.startswith(">"):
                if any_sequence:
                    yield determine_sequence_class(seq_header, any_sequence)
                    any_sequence = ""
                seq_header = line.strip().replace(">", "")
            else:
                any_sequence += line.strip()
        if any_sequence:
            yield determine_sequence_class(seq_header, any_sequence)

def determine_sequence_class(seq_header, sequence):
    """Determines if a sequence is DNA, RNA, or Protein and returns the appropriate instance."""
    try:
        if all(base in sequence_dictionaries.dna_letters for base in sequence):
            return DNASequence(seq_header, sequence)
        elif all(base in sequence_dictionaries.rna_letters for base in sequence):
            return RNASequence(seq_header, sequence)
        elif all(base in sequence_dictionaries.protein_letters for base in sequence):
            return ProteinSequence(seq_header, sequence)
        else:
            raise IncorrectSequenceLetter("Unknown letter", Sequence)
    except IncorrectSequenceLetter as e:
        sys.stderr.write(f"Skipping sequence {seq_header}: {e}\n")
        return None

class Sequence:
    alphabet = ""
    weights = {}

    def __init__(self, identifier, sequence):
        self.__identifier = identifier
        self.__sequence = sequence

        for residue in self.__sequence:
            if residue not in self.alphabet:
                raise IncorrectSequenceLetter(residue, self.__class__)

    def get_identifier(self):
        return self.__identifier

    def get_sequence(self):
        return self.__sequence

    def get_mw(self):
        return round(sum(self.weights.get(residue, 0) for residue in self.__sequence) - (len(self.__sequence) - 1) * 18.015, 2)

    def __len__(self):
        return len(self.get_sequence())

    def __str__(self):
        return f"{self.get_identifier()} == {self.get_sequence()}"

    def __lt__(self, other):
        return self.get_mw() < other.get_mw()

    def __hash__(self):
        return hash((self.get_identifier(), self.get_sequence()))

    def __eq__(self, other):
        return self.get_identifier() == other.get_identifier() and self.get_sequence() == other.get_sequence()

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
                translated_sequence += "*"  # Stop codon
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

def list_fasta_files(directory):
    """Lists all .fasta and .fa files in the specified directory."""
    fasta_files = []
    for filename in os.listdir(directory):
        if filename.endswith(".fasta") or filename.endswith(".fa"):
            fasta_files.append(os.path.join(directory, filename))
    return fasta_files

def process_files(input_path, output_file=None):
    """Processes input files and handles DNA, RNA, and Protein sequences."""
    if os.path.isdir(input_path):
        fasta_files = list_fasta_files(input_path)
    else:
        fasta_files = [input_path] if input_path.endswith((".fasta", ".fa")) else []

    sequences = []
    sys.stderr.write(f"{len(fasta_files)} FASTA files found.\n")

    for fasta_file in fasta_files:
        sys.stderr.write(f"Processing {fasta_file}...\n")
        for sequence in FASTA_iterator(fasta_file):
            if sequence:
                try:
                    if isinstance(sequence, (DNASequence, RNASequence)):
                        sequence = sequence.translate()  # Translate DNA/RNA to ProteinSequence
                    sequences.append(sequence)
                except IncorrectSequenceLetter as e:
                    sys.stderr.write(f"Skipping sequence {sequence.get_identifier()}: {e}\n")
        sys.stderr.write(f"{fasta_file} finished.\n")

    sys.stderr.write(f"{len(sequences)} sequences found.\n")
    sys.stderr.write("Sorting the sequences...Thanks for your patience.\n")
    sequences.sort()
    sys.stderr.write("Sort process finished.\n")

    output = "\n".join(f"{seq.get_identifier()}\t{len(seq)}\t{seq.get_mw()}" for seq in sequences)

    if output_file:
        with open(output_file, "w") as f:
            f.write(output + "\n")
    else:
        print(output)

    sys.stderr.write("Program finished correctly.\n")

if __name__ == "__main__":
    if len(sys.argv) == 1:
        process_files(".")
    elif len(sys.argv) == 2:
        process_files(sys.argv[1])
    elif len(sys.argv) == 3:
        process_files(sys.argv[1], sys.argv[2])
    else:
        sys.stderr.write("Use like: python python_exercise_block2_part4.py [IN] [OUT]\n")





