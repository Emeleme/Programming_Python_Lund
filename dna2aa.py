#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""

Name:dna2aa.py
Author:Maria Laura Mahecha Escobar
Description: Take a DNA fasta file and convert it into an aminoacid sequence
Institution:Lund University-Programming in python course
email:ma5668ma-s@student.lu.se
Date:2023y10m12d
Version:Python3.0 or more

User defined functions
- 

Procedure:
    1. Get DNA sequence in Fasta format.
    2. Check for errors in the input file
    3. Change T-Us in sequences in input file
    4. Read 3 bases to aminoacid code
    5. Save the new output file

Input: Fasta format DNA sequences
Output: Fasta format aminoacid sequences
    
Usage: 

    python dna2aa.py <DNA.fna> <output_file_name.txt>[optional]


"""

#####libraries#####

import os 
import sys


###Basic control check###
print("Checking input files: ...")

#for input length
####ESTO SIGUE SIN FUNCIONAR
# =============================================================================
# if 2 < len(sys.argv) > 3:
#     raise IndexError()
#     print("Please run this script in the form: python dna2aa.py <DNA.fna> <output_file_name.txt>[optional]")
#     sys.exit(1) #when other, then exit
# =============================================================================

if len(sys.argv) == 3: #when user inputs the output file
    output_file = open(sys.argv[2], "w") #output file is the user's name    
    print("The output of this process is going to be stored the file named: ", sys.argv[2])
elif len(sys.argv) == 2: #when user inputs just one file
    output_file = open("aminoacid_sequences.txt","w") #output file is aminoacid_sequences.txt
    print("The output of this process is going to be stored the file named: aminoacid_sequences.txt")

fasta_user=sys.argv[1]
print("The fasta file that you provided is: ", sys.argv[1])

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
        valid_nucleotides = ["A", "C", "G", "T"] # List of valid nucleotides
        sequence=[]
        for line in fasta_file: #for each line in the fasta file
            if not line.startswith(">"): #for each line that is not a header
                sequence += line.rstrip()
                sequence=str(line) #save the sequence in sequence variable as str
                sequence = sequence.upper().strip() #deletes spaces and \n chr
                for nucleotide in sequence: #for each nucleotide in the sequence
                    if nucleotide not in valid_nucleotides: #check if it is not in the valid nucleotides variable
                        raise TypeError (f"Seems that the fasta file you uploaded contain something different from nitrogenated bases. Please check your file and provide a file with only DNA sequences and no spaces at the end of the sequences ' '. Your file contains a '{nucleotide}'")
                        sys.exit(1) #if not raise an error and exit
        print("DNA Fasta file seems to be correct")

# =============================================================================
# for line in r:
#     if line.startswith('>'):
#         #Found start of a sequence
#         if temp:
#             seq.append(temp); #Keep the previous sequence in seq
#         temp='' #zero temp
#     else:
#         #Found more of existing sequence
#         temp += line.rstrip() #remove new line character
# 
# if temp:
#     #if the file is not empty
#     seq.append(temp)
# 
# r.close()
# =============================================================================

#For output file

#HACER!!! SI NO HAY UN NOMBRE HAGA UNO DEFAULT QUE SEA AMINOACIDS.TXT O ALGO ASI


#############SI EL FASTA FILE ESTÁ EN MÁS DE 1 LINEA VOLVERLO SOLO UNA LINEA
                
#####Commands#####

#Define the aminoacids dictionary for each codon. 
codontab = {
    'UCA': 'S',    # Ser
    'UCC': 'S',    # Ser
    'UCG': 'S',    # Ser
    'UCU': 'S',    # Ser
    'UUC': 'F',    # Phe
    'UUU': 'F',    # Phe
    'UUA': 'L',    # Leu
    'UUG': 'L',    # Leu
    'UAC': 'Y',    # Tyr
    'UAU': 'Y',    # Tyr
    'UAA': '*',    # Stop
    'UAG': '*',    # Stop
    'UGC': 'C',    # Cys
    'UGU': 'C',    # Cys
    'UGA': '*',    # Stop
    'UGG': 'W',    # Trp
    'CUA': 'L',    # Leu
    'CUC': 'L',    # Leu
    'CUG': 'L',    # Leu
    'CUU': 'L',    # Leu
    'CCA': 'P',    # Pro
    'CCC': 'P',    # Pro
    'CCG': 'P',    # Pro
    'CCU': 'P',    # Pro
    'CAC': 'H',    # His
    'CAU': 'H',    # His
    'CAA': 'Q',    # Gln
    'CAG': 'Q',    # Gln
    'CGA': 'R',    # Arg
    'CGC': 'R',    # Arg
    'CGG': 'R',    # Arg
    'CGU': 'R',    # Arg
    'AUA': 'I',    # Ile
    'AUC': 'I',    # Ile
    'AUU': 'I',    # Ile
    'AUG': 'M',    # Met
    'ACA': 'T',    # Thr
    'ACC': 'T',    # Thr
    'ACG': 'T',    # Thr
    'ACU': 'T',    # Thr
    'AAC': 'N',    # Asn
    'AAU': 'N',    # Asn
    'AAA': 'K',    # Lys
    'AAG': 'K',    # Lys
    'AGC': 'S',    # Ser
    'AGU': 'S',    # Ser
    'AGA': 'R',    # Arg
    'AGG': 'R',    # Arg
    'GUA': 'V',    # Val
    'GUC': 'V',    # Val
    'GUG': 'V',    # Val
    'GUU': 'V',    # Val
    'GCA': 'A',    # Ala
    'GCC': 'A',    # Ala
    'GCG': 'A',    # Ala
    'GCU': 'A',    # Ala
    'GAC': 'D',    # Asp
    'GAU': 'D',    # Asp
    'GAA': 'E',    # Glu
    'GAG': 'E',    # Glu
    'GGA': 'G',    # Gly
    'GGC': 'G',    # Gly
    'GGG': 'G',    # Gly
    'GGU': 'G'     # Gly
}
#codontab["GGA"]


new_fasta_file = output_file

print("Guetting amino acid sequences ...")

#To get aminoacid sequences from a DNA strand
with open (fasta_user, "r") as fasta_file: #open fasta file
    temp_sequences=fasta_file.readlines() #place all contents of fasta file in temp_sequences
    for line in temp_sequences: #for each line in temp_sequences
        if line.startswith (">"): #if the line is a header
            new_line = "header here" #changes the status of new_line
            line = line.strip() #takes the new line chr out
            new_fasta_file.write (line + "\n") #writes the header in the output file
        elif new_line is not None: #To add the next line in the output file, check if there was something in the new_line variable
            rna_seq = line.upper().strip() #store in a variable the sequence line uppercase and without spaces
            rna_seq = rna_seq.replace("T","U") #replace all T's in the sequence for U's and store it 
            aa="" #make a variable to add aminoacids (aa) on the go
            for i in range(0, len(rna_seq), 3): #for each third position in a range from 0 to the len of rna sequence, 
                cod = codontab[rna_seq[i:i+3]] #write a cod variable with the dictionary value for a 3 letters key starting with i
                aa += cod #add the value of the dictinary (aminoacid letter) to the aa variable
            new_fasta_file.write (aa + "\n") #write the complete aminoacid sequence in the output file
            new_line=None #set to None again to continue with the next line
            
print("Process complete")
print("Your output file is ready.")

#Close all files
fasta_file.close()
output_file.close()
