#########   Assignment_7: Object-Oriented Protein Analysis   ############
## * This script defines a Protein class to encapsulate protein data.      *
## * It includes a FASTA generator that yields Protein objects instead     *
## * of simple strings, allowing for complex biological calculations       *
## * like molecular weight and subsequence searching.                      *

import sys
# Ensure the file sequence_dictionaries.py is in the same folder
from sequence_dictionaries import protein_weights

class Protein():
    """Represents a biological protein sequence with analysis methods."""
    
    def __init__(self, identifier, sequence):
        # Attributes: The data the object holds
        self.identifier = identifier
        self.sequence = sequence
    
    def get_identifier(self):
        """Returns the protein ID."""
        return self.identifier
    
    def get_sequence(self):
        """Returns the amino acid sequence."""
        return self.sequence 
    
    def get_mw(self): 
        """
        Calculates Molecular Weight (Da).
        Includes water loss (18.015) for each peptide bond.
        """
        # Sum individual amino acids
        aa_sum = sum([protein_weights.get(residue, 0) for residue in self.sequence])
        # Subtract water for (n-1) bonds
        if len(self.sequence) > 1:
            total_mw = aa_sum - (len(self.sequence) - 1) * 18.015
        else:
            total_mw = aa_sum
        return round(total_mw, 2)
    
    def has_subsequence(self, other_protein):
        """Checks if another Protein object's sequence is inside this one."""
        return other_protein.get_sequence() in self.sequence
    
    def get_length(self):
        """Returns number of amino acids."""
        return len(self.sequence)

def FASTA_iterator(fasta_filename):
    """
    Generator: Reads FASTA and yields instances of the Protein class.
    Corrects the variable shadowing bug found in the previous version.
    """
    with open(fasta_filename , "r") as prot_file:
        current_sequence = ""
        current_id = None
        
        for line in prot_file:
            line = line.strip()
            if not line: continue # Skip empty lines

            if line.startswith(">"):                                                  
                # If we already have a stored sequence, yield the previous protein
                if current_id and current_sequence:
                    yield Protein(identifier=current_id, sequence=current_sequence)
                
                # Reset for the new protein record
                current_id = line[1:] # Remove '>' and save ID
                current_sequence = ""                                                 
            else:                                                                       
                current_sequence += line                                              
        
        # Yield the very last protein in the file
        if current_id and current_sequence:
            yield Protein(identifier=current_id, sequence=current_sequence)

# --- Execution Block ---
if __name__ == "__main__":
    # Check if a filename was provided in terminal
    if len(sys.argv) > 1:
        fasta_input = sys.argv[1]
        
        # We iterate through the generator
        for protein_obj in FASTA_iterator(fasta_input):
            # Print in the format requested by the professor
            print(protein_obj.get_length())
            print("%s (%d): %.4f" % (
                protein_obj.get_identifier(), 
                protein_obj.get_length(), 
                protein_obj.get_mw()
            ))
    else:
        print("Usage: python script.py <fasta_file>")


#########   Assignment_7: Object-Oriented Protein Analysis   ############ RESULTS
##Call### python3 python_exercise_block2_part1.py example_fasta_file.fa 

##Output#
#1732 
# Q5VT25 (1732): 197308.1500
#
#1732: This is the length of the protein sequence in amino acids. The "Raw" length.
#Q5VT25 (Identifier): This is the UniProt Accession number. It functions like a unique ID card for a specific protein.
#(1732) (Length): This tells you the protein consists of 1,732 amino acids.
#197308.1500 (Molecular Weight): This is the mass of the protein in Daltons (Da). In biochemistry, we often convert this to kiloDaltons (kDa) for readability, so this protein is approximately 197.3 kDa.