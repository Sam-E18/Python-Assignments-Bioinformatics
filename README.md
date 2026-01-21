# Python-Assignments-Bioinformatics 🧬

This repository contains assignments from the **Introduction to Python** course, specifically tailored for **Bioinformatics in Health Sciences**. Most scripts focus on the biological analysis of sequences, molecular weights, and data parsing.

> [!IMPORTANT]
> **Note:** These exercises are intended as a guide. Some scripts include personal modifications or scientific improvements (like water-loss calculations in proteins) that may differ from the basic assignment requirements. Always verify the code logic before use.

---

## 📂 Summary of Assignments

* **Assignment 1 - 6:** Foundational Python logic, including string manipulation, list comprehension, and basic dictionary usage for biological data.
* **Assignment 7 (Object-Oriented Proteins):** * Defines a `Protein` class to encapsulate identifiers and sequences.
    * Implements a custom **FASTA Generator** to handle large genomic files.
    * Calculates scientifically accurate **Molecular Weight** (accounting for water loss).
* **Assignment 8 (Class Hierarchy & Inheritance):** * Creates a `Sequence` superclass with specialized subclasses: `DNASequence`, `RNASequence`, and `ProteinSequence`.
    * Automates biological workflows: **DNA ⮕ RNA (Transcription)** and **RNA ⮕ Protein (Translation)**.
    * Uses polymorphism to manage different molecular weight calculations across types.

---

## 🚀 Workflow Suggestions for Beginners

I am not an expert, but here are the suggestions I wish I had followed when I started:

* **Editor:** Download **Visual Studio Code (VSC)**. It is extremely beginner-friendly and has a massive library of extensions.
    * *Tip:* Use the "vscode-pdf" extension to read assignment instructions side-by-side with your code.
* **Environment:** Install Python locally and link it to VSC. This makes managing scripts and libraries (like `sys` or `math`) much easier.
* **Structure:** Keep a `sequence_dictionaries.py` file in your root folder. It acts as a "central library" for your DNA tables and amino acid weights.
* **Testing:** Always test your scripts using the `if __name__ == "__main__":` block to ensure they run correctly from the terminal.

---

## 🛠️ Technical Stack
* **Language:** Python 3.x
* **Libraries used:** `sys` (for terminal arguments), `math`, and custom modules.
* **Key Concepts:** Object-Oriented Programming (OOP), Generators, Class Inheritance, and Data Parsing.

---
*Created with ❤️ for the Bioinformatics community.*