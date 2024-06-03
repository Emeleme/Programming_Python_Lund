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
- one_line_seq(x): concatenate all the sequences of a fasta file in one line

Procedure:
    1. Read a fasta file containing various sequences and make it just one line
sequence
    2. Go through all the sequence and count each aminoacid abundance
    3. Print the global abundances of the aminoacids in a new file

Input: Fasta format protein sequences
Output: .txt file containing complete aminoacids counts
    
Usage: 

    python amino_count.py <protein.fna> <output_file_name.txt>[optional]

"""

#####libraries#####

import os 
import sys


###Basic control check###
print("Checking input files: ...")

#for input length
if len(sys.argv) == 3: #when user inputs the output file
    fasta_user=sys.argv[1] #first argument is fasta file
    print("The fasta file that you provided is: ", sys.argv[1])
    output_file = open(sys.argv[2], "w") #output file is the user's name    
    print("The output of this process is going to be stored the file named: ", sys.argv[2])
elif len(sys.argv) == 2: #when user inputs just one file
    fasta_user=sys.argv[1] #first argument is fasta file
    print("The fasta file that you provided is: ", sys.argv[1])
    output_file = open("aminoacid_count.txt","w") #output file is aminoacid_sequences.txt
    print("The output of this process is going to be stored the file named: aminoacid_sequences.txt")
elif len(sys.argv) == 1: #when user just types the name of the code
    print("Please provide a fasta file and use this code like this: python dna2aa.py <DNA.fna> <output_file_name.txt>[optional]")
    sys.exit(1) #if not then exit
elif len(sys.argv) >= 4: #when the user types more than expected
    print("There are too many files, did you typed something extra? Please use this code like this: python dna2aa.py <DNA.fna> <output_file_name.txt>[optional]")
    sys.exit(1) #if not then exit
    
#for fasta file
if not os.path.exists(fasta_user): #Does the file exist?
    print("The fasta file you uploaded does not exist. Please provide an existing file.")
    sys.exit(1) #if not then exit
elif not os.path.getsize(fasta_user)>0: #Does the file contain values?
    print("The fasta file you uploades does not contain any sequence. Please provide a sequence")
    sys.exit(1) #if not then exit

with open (fasta_user, "r") as fasta_file:
    #chech if the file is a fasta file
    first_line_fasta=fasta_file.readline() #Does the first line of the fasta file starts with something different of ">"
    if not first_line_fasta.startswith (">"):
        raise TypeError("It seems you uploaded a non supported type of file, maybe a Fastq? Please check and provide a fasta format file")
        sys.exit(1) #if yes then exit
    else:
        #checking if the fasta file contains ONLY valid nucleotides
        valid_aminoacids = ["A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N", "O", "P", "Q", "R", "S", "T", "U", "V", "W", "X", "Y", "Z"] # List of valid letters
        sequence=[]
        for line in fasta_file: #for each line in the fasta file
            if not line.startswith(">"): #for each line that is not a header
                sequence += line.rstrip()
                sequence=str(line) #save the sequence in sequence variable as str
                sequence = sequence.upper().strip() #deletes spaces and \n chr
                for aminoacid in sequence: #for each aminoacid in the sequence
                    if aminoacid not in valid_aminoacids: #check if it is not in the valid nucleotides variable
                        raise TypeError (f"Seems that the fasta file you uploaded contain something different from letters. Please check your file and provide a file with only aminoacid sequences. Your file contains a '{aminoacid}'")
                        sys.exit(1) #if not raise an error and exit
        print("DNA Fasta file seems to be correct")





#Concatenate all the sequences in one line

#Make a function to store the sequences of one file in just one line
def one_line_seq(x):
    
    fasta_file = open (x, "r")
    file = fasta_file.readlines() #read all the file into a new variable
    
    sequence="" #open an empty variable
    for line in file: #for each line in the file variable
        if not line.startswith(">"): #if its not a header
            sequence += line.strip() #add all sequences without the new line chr
            sequence = sequence.upper() #make sequence uppercase
    fasta_file.close() #close file
    return sequence 

sequence = one_line_seq(fasta_user) #run function on fasta file from user   
    
aminoacids = ["A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y"] #Define aminoacid list
matching="" #open an empty variable for matching aminoacids
non_matching="" #open an empty variable for nonmatching aminoacids
for i in sequence: #for each position in the sequence
    if i in aminoacids: #if matches with the valid aminoacids
        letter = i #store the letter
        matching += letter #add the letter in matching variable
    else: #if does not metch
        other_letter = i #store the letter
        non_matching += other_letter #add the letter to non_matching variable


counts=[] #open an empty list
for j in aminoacids: #foe each letter in aminoacids variable
    counts += [j + " " + str(matching.count(j))] #count how many of that letter there is in matching variable and store the result as "letter count" in count variable

counts += ["X "+str(len(non_matching))] #add at the end of count variable the length of the nonmatching variable as "X"
  
print("Writing output file ...")  
aminoacid_output = output_file #set output file to new variable
for k in counts: #for item in counts
    aminoacid_output.write (k + "\n") #write in the output file the item + new line chr

print("Process complete")
#close file
aminoacid_output.close()
