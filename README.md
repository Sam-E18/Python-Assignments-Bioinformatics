# Python-Assignments-Bioinformatics 🧬

This repository contains assignments from the **Introduction to Python** course, specifically tailored for **Bioinformatics in Health Sciences**. Most scripts focus on the biological analysis of sequences, molecular weights, and data parsing.

> [!IMPORTANT]
> **Note:** These exercises are intended as a guide. Some scripts include personal modifications or scientific improvements (like water-loss calculations in proteins) that may differ from the basic assignment requirements. Always verify the code logic before use.

---

##  Summary of Assignments
 **Foundational Python logic, including string manipulation, list comprehension, and basic dictionary usage for biological data.**

* **Assignment 1: Environment Setup**
  * Introduction to the Python interpreter and basic variable assignments.
* **Assignment 2: String Manipulation**
   * Basic operations on DNA strings: counting nucleotides and finding simple motifs.
* **Assignment 3: Lists & GC Content**
   * Using lists to store sequence data and calculating the GC percentage of a DNA strand.
* **Assignment 4: Dictionaries & Codons**
   * Implementing dictionaries to map amino acids to their molecular weights and codons to residues.
* **Assignment 5: Basic File Reading**
   * Opening and reading simple text files containing single biological sequences.
* **Assignment 6: Advanced FASTA Parsing**
   * Developing logic to handle multi-line FASTA files and extracting headers from sequences.
* **Assignment 7 (Object-Oriented Proteins):** * Defines a `Protein` class to encapsulate identifiers and sequences.
    * Implements a custom **FASTA Generator** to handle large genomic files.
    * Calculates scientifically accurate **Molecular Weight** (accounting for water loss).
* **Assignment 8 (Class Hierarchy & Inheritance):** * Creates a `Sequence` superclass with specialized subclasses: `DNASequence`, `RNASequence`, and `ProteinSequence`.
    * Automates biological workflows: **DNA ⮕ RNA (Transcription)** and **RNA ⮕ Protein (Translation)**.
    * Uses polymorphism to manage different molecular weight calculations across types.
* **Assignment 9: Special Methods & Hashing** 
    * Implements "Special" methods (__len__, __add__, __str__, etc.) to make sequences behave like native Python objects.
    * Enables complex sorting based on Molecular Weight using comparison operators.
    * Implements __hash__ to allow Sequence objects to be used as keys in dictionaries or stored in sets.
* **Assignment 10: The Integrated Pipeline (Professional Tool)**
    * **Advanced CLI:** Implements argparse for industry-standard command-line flags (-i, -o, -v).

    * **Biological Precision:** Unlike previous versions, this translation logic mimics real biology by searching for Start Codons and terminating at Stop Codons.

    * **Performance Optimization:** Uses Python sets for O(1) residue validation and Caching to store Molecular Weight results, preventing redundant calculations.

    * **Robust Error Handling:** Features a custom exception class IncorrectSequenceLetter to skip corrupt data without crashing the entire pipeline.
---

##  Workflow Suggestions for Beginners

I am not an expert, but here are the suggestions I wish I had followed when I started:

* **Editor:** Download **Visual Studio Code (VSC)**. It is extremely beginner-friendly and has a massive library of extensions.
    * *Tip:* Use the "vscode-pdf" extension to read assignment instructions side-by-side with your code.
* **Environment:** Install Python locally and link it to VSC. This makes managing scripts and libraries (like `sys` or `math`) much easier.
* **Structure:** Keep a `sequence_dictionaries.py` file in your root folder. It acts as a "central library" for your DNA tables and amino acid weights.
* **Testing:** Always test your scripts using the `if __name__ == "__main__":` block to ensure they run correctly from the terminal.

---

##  Technical Stack
* **Language:** Python 3.x
* **Libraries used:** `sys` (for terminal arguments), `math`, and custom modules.
* **Key Concepts:** Object-Oriented Programming (OOP), Generators, Class Inheritance, and Data Parsing.

---
*Created with ❤️ for the Bioinformatics community.*
