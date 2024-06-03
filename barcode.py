#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""

Name:barcode.py
Author:Maria Laura Mahecha Escobar
Description: Save sequences with definded barcode in a new fastq file
Institution:Lund University-Programming in python course
email:ma5668ma-s@student.lu.se
Date:2023y10m15d
Version:Python3.0 or more

User defined functions
- 

Procedure:
    1. Read a fastq file sequences  
    2. Go through all the file reading sequences and saving the structure for each header
    3. Match end or the begining of the sequences with the barcodes
    4. Write matching sequences into one file
    5. Write non matching sequences into another file
    
Input: Fastq format DNA sequences
Output: two fastq files containing barcode sequences (sample1.fastq) and nonmatching sequences to barcodes (sample2.fastq)
    
Usage: 

    python barcode.py <DNA.fastq> <output_file_name1.txt>[optional] <output_file_name2.txt>[optional]

"""

#####libraries#####

import os 
import sys


###Basic control check###
print("Checking input files: ...")

#for input length
if len(sys.argv) == 4: #when user inputs the output file
    fastq_user=sys.argv[1] #first argument is fasta file
    print("The fastq file that you provided is: ", sys.argv[1])
    output_file1 = open(sys.argv[2], "w") #first output file is the user's name    
    print("The sequences finishing with specified barcodes are going to be stored in the file named: ", sys.argv[2])
    output_file2 = open(sys.argv[3] , "w") #second output file is sample_2.txt
    print("The sequences without barcode are going to be stored in the file named: ", sys.argv[3])
elif len(sys.argv) == 3: #when user inputs the output file
    fastq_user=sys.argv[1] #first argument is fasta file
    print("The fastq file that you provided is: ", sys.argv[1])
    output_file1 = open(sys.argv[2], "w") #first output file is the user's name    
    print("The sequences finishing with specified barcodes are going to be stored in the file named: ", sys.argv[2])
    output_file2 = open("sample2.fastq" , "w") #second output file is sample_2.txt
    print("The sequences without barcode are going to be stored in the file named: sample2.fastq")
elif len(sys.argv) == 2: #when user inputs just one file
    fastq_user=sys.argv[1] #first argument is fasta file
    print("The fastq file that you provided is: ", sys.argv[1])
    output_file1 = open("sample1.fastq", "w") #first output file is sample_1.txt    
    print("The sequences finishing with specified barcodes are going to be stored in the file named: sample1.fastq")
    output_file2 = open("sample2.fastq" , "w") #second output file is sample_2.txt
    print("The sequences without barcode are going to be stored in the file named: sample2.fastq")
elif len(sys.argv) == 1: #when user just types the name of the code
    print("Please provide a fastq file and use this code like this: python barcode.py <DNA.fastq> <output_file_name1.txt>[optional] <output_file_name2.txt>[optional]")
    sys.exit(1) #if not then exit
elif len(sys.argv) >= 5: #when the user types more than expected
    print("There are too many files, did you typed something extra? Please use this code like this: python barcode.py <DNA.fastq> <output_file_name1.txt>[optional] <output_file_name2.txt>[optional]")
    sys.exit(1) #if not then exit

#for fastq file
if not os.path.exists(fastq_user): #Does the file exist?
    print("The fastq file you uploaded does not exist. Please provide an existing file.")
    sys.exit(1) #if not then exit
elif not os.path.getsize(fastq_user)>0: #Does the file contain values?
    print("The fastq file you uploaded does not contain any sequence. Please provide a sequence")
    sys.exit(1) #if not then exit

with open (fastq_user, "r") as fastq_file:
    #chech if the file is a fasta file
    first_line_fastq=fastq_file.readline() #Does the first line of the fasta file starts with something different of ">"
    if not first_line_fastq.startswith ("@"):
        raise TypeError("It seems you uploaded a non supported type of file, maybe a Fasta file? Please check and provide a fastq format file")
        sys.exit(1) #if yes then exit
    else:
        print("Fastq file seems to be correct")


#####Commands#####
fastq_file = open(fastq_user, "r") #open fastq file of the user

print("writing output files ...")
for line in fastq_file: #for echa line in the fastq file
    header1 = line.rstrip() #save the line (header) in a variable without new line chr
    sequences = next(fastq_file).rstrip() #save the next line (sequence line) in another variable without new line chr
    header2 = next(fastq_file).rstrip() #save the next line (+ line) in a variable called header2 without new line chr
    qualities = next(fastq_file).rstrip() #save the next line (quality scores line) in a variable without new line chr
    if sequences.startswith("TATCCTCT" or "GTAAGGAG" or "TCTCTCCG"): #if the beggining of the sequence is any barcode
        sequences = sequences[:-8] #remove the barcode from the sequence
        qualities = qualities[:-8] #remove the quality scores for the sequence as well
        output_file1.write(header1 + "\n" + sequences + "\n" + header2 + "\n" + qualities + "\n") #write the whole sequence structure in the sample1 file
    elif sequences.endswith("TATCCTCT" or "GTAAGGAG" or "TCTCTCCG"): #if the end of the sequence is any barcode
        sequences = sequences[9:] #remove the barcode from the sequence
        qualities = qualities[9:] #remove the quality scores for the sequence as well
        output_file1.write(header1 + "\n" + sequences + "\n" + header2 + "\n" + qualities + "\n") #write the whole sequence structure in the sample1 file
    else:
        output_file2.write(header1 + "\n" + sequences + "\n" + header2 + "\n" + qualities + "\n") #if the sequence does not have a barcode nor at the begining or the end, write the structure in the sample2 file
  
print("Process complete")

#close all files
fastq_file.close()
output_file1.close()
output_file2.close()

