#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""

Name:amino_count.py
Author:Maria Laura Mahecha Escobar
Description: Count amino acid abundances in a file
Institution:Lund University-Programming in python course
email:ma5668ma-s@student.lu.se
Date:2023y10m13d
Version:Python3.0 or more

User defined functions
- 

Procedure:
    1. Read a fasta file containing various sequences and make it just one line
sequence
    2. Go through all the sequence and count each aminoacid 
    3. Print the results

Input: Fasta format protein sequences
Output: .txt file containing individual aminoacids counts
    
Usage: 

    python amino_count.py <protein.fna> <output_file_name.txt>[optional]


"""

#CHECK THAT THE AMINOACID SEQUENCE HAS ONLY LETTERS
#FALTA HACER TODO EL CHEQUEO DE CALIDAD

#Concatenate all the sequences in one line

#Make a function to store the sequences of one file in just one line
def one_line_seq(x):
    
    fasta_file = open (x, "r")
    file = fasta_file.readlines()
    
    sequence=""
    for line in file:
        if not line.startswith(">"):
            sequence += line.strip()
            sequence = sequence.upper()
    fasta_file.close()
    return sequence

sequence = one_line_seq("aminoacids_seq.txt")       
    
aminoacids = ["A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y"]
matching=""
non_matching=""
for i in sequence:
    if i in aminoacids:
        letter = i
        matching += letter
    else:
        other_letter = i
        non_matching += other_letter


counts=[]
for j in aminoacids:
    counts += [j + " " + str(matching.count(j))]

counts += ["X "+str(len(non_matching))]
    
    
with open ("aminoacid_count.txt", "w") as aminoacid_output:
    for k in counts:
        aminoacid_output.write (k + "\n")

#close file
aminoacid_output.close()
