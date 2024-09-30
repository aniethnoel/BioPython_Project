# BioPython_Project
# Introduction
This project aims to analyze the gene with accession number AB000824.1 using Biopython and other bioinformatics tools. The tasks involve gene sequence retrieval, transcription, translation, primer design, domain extraction, homolog analysis, and multiple sequence alignment (MSA). The key components are implemented in Python, leveraging the Biopython library.

# Features
1. Download Gene Sequence: Download the gene sequence associated with the accession number and save it in FASTA format.
2. mRNA Conversion: Convert the downloaded gene sequence to mRNA.
3. Protein Translation: Translate the mRNA sequence into a protein sequence, varying the open reading frame (ORF).
4. Primer Design: Design forward and reverse primers of 30 nucleotides for the gene.
5. Protein Sequence Retrieval: Retrieve the protein sequence from the GenBank file based on the protein ID.
6. Catalytic Domain Extraction: Extract the catalytic domain from the protein sequence using InterProScan for domain information.
7. BLASTP for Domain Analysis: Perform a BLASTP search with the extracted domain sequence.
8. Homologous Sequences: Retrieve five homologous sequences from different organisms, change their accession IDs to organism names.
9. Multiple Sequence Alignment (MSA): Perform MSA using Clustal Omega and interpret the relationships between organisms based on sequence similarity.

# Requirements
**Python Packages:**
1. Biopython
2. Requests (for sequence download)
3. Clustal Omega (for MSA)
5. BLAST+ (for BLASTP analysis)

**External Tools:**
1. InterProScan (for domain extraction)
2. BLASTP (for homolog sequence search)

# Results Interpretation
After running the MSA, examine the resulting aligned sequences and identify closely related organisms based on sequence similarity.
