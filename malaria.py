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
fasta_file=open("malaria.fna", "r")
blast_file=open("malaria.blastx.tab", "r")
new_fasta_file=open("output.txt","w")

#reading files
line1=fasta_file.readline()
print(line1)
print(line1.split())
for line in fasta_file:
    print(line)
    

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