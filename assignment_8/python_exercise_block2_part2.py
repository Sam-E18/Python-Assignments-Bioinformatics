######### Assignment_8: Sequence Class Hierarchy ############
## * This script implements a class hierarchy for biological sequences. *
## * It uses inheritance to define specific behaviors for DNA, RNA, and *
## * Proteins, allowing for operations like transcription and translation. *
#######################################################################

import sequence_dictionaries 

class Sequence():
    "Super class Sequence to define a sequence, there are sub clases into the super class"
    alphabet = ""
    weights = {}
    def __init__(self, identifier, sequence ): ## init the sequence instance 
        self.__identifier = identifier
        self.__sequence = sequence

    def get_identifier(self): ## method to obtain the identifier 
        return self.__identifier
    
    def get_sequence(self): ## method to get the sequence string
        return self.__sequence 

    
    def get_mw(self): ## its better to put the calculation in one line with the generator.
        ## method to obtain the meanweigth
        return round(sum(self.weights.get(residue, 0) for residue in self.__sequence) - (len(self.__sequence) - 1) * 18.015, 2)

    
    def has_subsequence(self, sequence_objt): ## this could be class Sequence or ProteinSequence or DNASequence or RNAsequence.
        return sequence_objt.get_sequence() in self.get_sequence()
    
    

class ProteinSequence(Sequence):
    """subclass to define a protein into the sequence"""
    def __init__(self, identifier, sequence):
        super().__init__(identifier, sequence)
    alphabet = sequence_dictionaries.protein_letters
    weights = sequence_dictionaries.protein_weights
   
    

class NucleotideSequence(Sequence):
    """ subclass nucleotidesequence to define nucletides into sequence"""
    def __init__(self, identifier, sequence):
        super().__init__(identifier, sequence)

    translate_table = {}

    def translate(self): 
        translated_seq = ""
        for codon_start in range(0,len(self.get_sequence()),3):
            codon = self.get_sequence()[codon_start:codon_start+3]
            translated_seq += self.translate_table[codon]
        translated_seq = ProteinSequence(self.get_identifier(), translated_seq)
        return translated_seq
    

class DNASequence(NucleotideSequence):
    """ subclass DNAsequence into Nucleotidesequence"""
    translate_table = sequence_dictionaries.dna_table
    alphabet = sequence_dictionaries.dna_letters
    weights = sequence_dictionaries.dna_weights

    def transcribe(self):
        transcribed_seq = ""
        for nucleobase in self.get_sequence():
            if nucleobase == "T":
                transcribed_seq += "U"
            else: 
                transcribed_seq += nucleobase
        transcribed_seq = RNASequence(self.get_identifier(), transcribed_seq)
        return transcribed_seq


class RNASequence(NucleotideSequence):
    """ subclass RNAsequence into Nucleotidesequence"""
    translate_table = sequence_dictionaries.rna_table 
    alphabet = sequence_dictionaries.rna_letters 
    weights = sequence_dictionaries.rna_weights

    def reverse_transcribe(self):
        retrieve_sequence = ""
        for nucleobase in self.get_sequence():
            if nucleobase == "U":
                retrieve_sequence += "T"
            else:
                retrieve_sequence += nucleobase
        retrieve_sequence = DNASequence(self.get_identifier() , retrieve_sequence)

        return retrieve_sequence
    

## define the instances only for Proteine, dna, and rna
p = ProteinSequence("protein1", "ACDEFGG")
dna = DNASequence("dna1", "ATGCATGC")
rna = RNASequence("rna1", "AUGCAUGC")

## print the instances 
print(p.get_identifier())  
print(p.get_sequence())   
print(p.get_mw())          #print the corrected weights

print(dna.get_sequence())  

print(dna.get_identifier())  
print(dna.get_sequence())    
print(dna.transcribe().get_sequence())  

print(rna.reverse_transcribe().get_sequence())  

print(p.has_subsequence(ProteinSequence("P2", "CDEF")))  
print(dna.has_subsequence(DNASequence("D2", "GCAT")))  




#####################################################################################

import sequence_dictionaries 

class Sequence():
    "Super class Sequence to define a sequence, there are sub clases into the super class"
    alphabet = ""
    weights = {}
    
    def __init__(self, identifier, sequence ): ## init the sequence instance 
        self.__identifier = identifier
        self.__sequence = sequence.upper() # Added .upper() to ensure logic works with dictionaries

    def get_identifier(self): ## method to obtain the identifier 
        return self.__identifier
    
    def get_sequence(self): ## method to get the sequence string
        return self.__sequence 

    def get_mw(self): ## its better to put the calculation in one line with the generator.
        ## method to obtain the meanweigth
        # Base calculation: just the sum of individual weights
        return sum(self.weights.get(residue, 0) for residue in self.__sequence)

    def has_subsequence(self, sequence_objt): ## this could be class Sequence or ProteinSequence or DNASequence or RNAsequence.
        return sequence_objt.get_sequence() in self.get_sequence()
    

class ProteinSequence(Sequence):
    """subclass to define a protein into the sequence"""
    def __init__(self, identifier, sequence):
        super().__init__(identifier, sequence)
    
    alphabet = sequence_dictionaries.protein_letters
    weights = sequence_dictionaries.protein_weights

    def get_mw(self): 
        """
        Overrides the base get_mw to include the water loss logic 
        specific to peptide bonds in proteins.
        """
        aa_weights = super().get_mw()
        # Proteins lose 18.015 Da for every bond (n-1)
        if len(self.get_sequence()) > 1:
            return round(aa_weights - (len(self.get_sequence()) - 1) * 18.015, 2)
        return round(aa_weights, 2)
   

class NucleotideSequence(Sequence):
    """ subclass nucleotidesequence to define nucletides into sequence"""
    def __init__(self, identifier, sequence):
        super().__init__(identifier, sequence)

    translate_table = {}

    def translate(self): 
        translated_seq = ""
        # Improved range to ensure we only process complete triplets (codons)
        for codon_start in range(0, len(self.get_sequence()) - (len(self.get_sequence()) % 3), 3):
            codon = self.get_sequence()[codon_start : codon_start + 3]
            translated_seq += self.translate_table.get(codon, "") # .get() prevents crashing on unknown codons
        
        # Return a new ProteinSequence object
        translated_seq = ProteinSequence(self.get_identifier(), translated_seq)
        return translated_seq
    

class DNASequence(NucleotideSequence):
    """ subclass DNAsequence into Nucleotidesequence"""
    translate_table = sequence_dictionaries.dna_table
    alphabet = sequence_dictionaries.dna_letters
    weights = sequence_dictionaries.dna_weights

    def transcribe(self):
        # Using .replace() is much faster than a loop, but keeps your variable name
        transcribed_seq = self.get_sequence().replace("T", "U")
        
        # Return a new RNASequence object
        transcribed_seq = RNASequence(self.get_identifier(), transcribed_seq)
        return transcribed_seq


class RNASequence(NucleotideSequence):
    """ subclass RNAsequence into Nucleotidesequence"""
    translate_table = sequence_dictionaries.rna_table 
    alphabet = sequence_dictionaries.rna_letters 
    weights = sequence_dictionaries.rna_weights

    def reverse_transcribe(self):
        # Using .replace() is more efficient for biological sequences
        retrieve_sequence = self.get_sequence().replace("U", "T")
        
        # Return a new DNASequence object
        retrieve_sequence = DNASequence(self.get_identifier() , retrieve_sequence)
        return retrieve_sequence
    

# --- Testing the Instances ---
if __name__ == "__main__":
    ## define the instances only for Proteine, dna, and rna
    p = ProteinSequence("protein1", "ACDEFGG")
    dna = DNASequence("dna1", "ATGCATGC")
    rna = RNASequence("rna1", "AUGCAUGC")

    ## print the instances 
    print(f"Protein ID: {p.get_identifier()}")  
    print(f"Protein Seq: {p.get_sequence()}")   
    print(f"Protein MW: {p.get_mw()}") # Correctly subtracts water

    print(f"DNA ID: {dna.get_identifier()}")  
    print(f"DNA Seq: {dna.get_sequence()}")    
    print(f"DNA MW: {dna.get_mw()}") # Does NOT subtract water (biologically correct)
    
    # Show Transcription
    transcribed = dna.transcribe()
    print(f"Transcribed RNA: {transcribed.get_sequence()}")  

    # Show Translation (from DNA directly)
    translated = dna.translate()
    print(f"Translated Protein: {translated.get_sequence()}")

    # Subsequence checks
    print(f"Protein has 'CDEF'?: {p.has_subsequence(ProteinSequence('P2', 'CDEF'))}")  
    print(f"DNA has 'GCAT'?: {dna.has_subsequence(DNASequence('D2', 'GCAT'))}")