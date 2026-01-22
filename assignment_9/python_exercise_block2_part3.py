import sequence_dictionaries 

class Sequence():
    """ Define the differents behaviour for the classes defined in the previous exercise block 2 part 2"""
    alphabet = ""
    weights = {}
    def __init__(self, identifier, sequence ):
        self.__identifier = identifier
        self.__sequence = sequence

        for residue_prot_dna_rna in self.__sequence: ##any of the residues form prot dna and rna
            if residue_prot_dna_rna not in self.alphabet:
                raise ValueError(f"Imposssible to create instance: {residue_prot_dna_rna} not possible") #show error if the instances is not dna rna or prot.
                
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
            raise TypeError("only one type of sequence") ## this is part of the other home work but its work better with that here.
        
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
    def __init__(self, identifier, sequence):
        super().__init__(identifier, sequence)
    alphabet = sequence_dictionaries.protein_letters
    weights = sequence_dictionaries.protein_weights
   
    

class NucleotideSequence(Sequence):
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
dna0 = DNASequence("dna0", "ATGCATGCTTTTT")
dna1 = DNASequence("dna1", "ATCG")
dna2 = DNASequence("dna2", "ATGC")
dna3 = dna0 + dna1 
p = ProteinSequence("P1", "ACDEFGG")
rna = RNASequence("R1", "AUGCAUGC")
dna_duplicate = DNASequence("dna0", "ATGCATGCTTTTT")  # Same identifier and sequence as dna0
# Use in a dictionary
sequence_dictionary = {
    dna0: "DNA 0 description",
    dna1: "DNA 1 description"
}
# Use in a set
sequence_set = {dna0, dna1, dna_duplicate}

print("Length of the sequence: %s" %len(dna0)) ## print the length
print(dna1 == dna2) ## should print false
print(dna1 != dna0) ## should print true
print(dna3) ## print concatenate sequence
print(dna0[0:3])
print(dna2 in dna0)  # Should return True because "ATGC" is in "ATGCATGCTTTTT"

print(f"dna1 < dna0: {dna1 < dna0}") ### __lt__ (less than), True if the molecular weight of dna1 is less than that of dna0.
print(f"dna1 <= dna0: {dna1 <= dna0}") ## __le__ (less than or equal to), True if the molecular weight of dna1 is less than or equal to that of dna0.
print(f"p > rna: {p > rna}")  ## __gt__ (greater than), True if the molecular weight of p is greater than that of rna.
print(f"p >= dna0: {p >= dna0}")  ## __ge__ (greater than or equal to), True if the molecular weight of p is greater than or equal to that of dna0.

## comparing sequence 
sequences = [dna0, dna1, p, rna]
sorted_sequences = sorted(sequences)
for seq in sorted_sequences:
    print(f"{seq.get_identifier()}: MW = {seq.get_mw()} g/mol")

#method hash
print(sequence_dictionary[dna0])               # Output: "DNA 0 description"
print(dna_duplicate in sequence_dictionary)    # Output: True, because dna_duplicate is considered the same as dna0
print(len(sequence_set))  # Output should be 2, because dna0 and dna_duplicate are considered the same
