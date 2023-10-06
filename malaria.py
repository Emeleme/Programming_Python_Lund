#!/usr/bin/env python
# coding: utf-8

###About###

#Name:malaria.py
#Author:Maria Laura Mahecha Escobar
#Description:Script to add the BLAST hit description to the fasta sequence header
#Institution:Lund University-Programming in python course
#email:ma5668ma-s@student.lu.se
#Date:2023y10m06d
#Version:Python3.0 or more

###Libraries### 
import sys


###Comands### 
#Opening files
fasta_file=open("five_seq.fna", "r")
blast_file=open("five_seq.blastx.tab", "r")
new_fasta_file=open("output.txt","w")


    
# =============================================================================
# #Check if the sequences have the right characters ALGO ESTÁ MAL CON ESTO. MIRAR MAS EXCEPCIONES 
# valid_nucleotides=["A","T","C","G"] #This are the only nucleotides accepted
# fasta_check=[] #open a new variable to check fasta file sequences
# for line in fasta_file:
#     if line[0] != ">":
#         print(line)
#         seq_check=line
#         fasta_check.append(seq_check)
# print(fasta_check)
# for nucleotide in str(fasta_check).upper(): #for each position in the sequence provided, check if it is in the valid values provided and then decide if theres an exception or not
#     if nucleotide not in valid_nucleotides: #If the letters in the 
#         raise Exception("Not valid nucleotides provided")
# else:
#     print("The file you uploaded contains DNA sequences")
# =============================================================================




#For each line in fasta_file search for lines begining with ">" and extract all
#information until another the next ">" appears and store the result in vasiable
#named id_seq 

seq_fas=[] #open a variable to store temporarly all the fasta sequences
new_fasta=[] #open a new variable to store all the modified sequences in fasta format

for line in fasta_file:
    if line.startswith(">"):
        if seq_fasta:
            seq.append(temp); #Keep the previous sequence in seq
        temp=''















#Closing all files
fasta_file.close()
blast_file.close()
new_fasta_file.close()