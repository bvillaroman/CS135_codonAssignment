# CS135_codonAssignment

This program was written for an assignment in a Software Analysis and Design class and was written without knowledge of how to use classes in C++.

=======================
 Reads a txt file called "dna.txt" containing invalid and valid DNA sequences and creates three new text files called:

##### headerLinesFreq
  - contains the frequency of the number of header lines in dna.txt

##### GCcontent
  - contains the original DNA data, sorted by gc content that lists gc content of each valid DNA sequence along with their header line

##### expandedDNA
  - contains everything in GCcontent
  - along with their transcription of DNA into RNA
  - the number of times the sequence "UGG" exists in the RNA sequence
  - the complement DNA string 
