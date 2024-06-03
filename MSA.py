#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Name:MSA.py
Author:Maria Laura Mahecha Escobar
Description: Reads the fasta file of multiple alignment sequences and returns 
their identity score and the overall score of the alignment
Institution:Lund University-Programming in python course
email:ma5668ma-s@student.lu.se
Date:2023y10m22d
Version:Python3.0 or more

User defined functions: None

Procedure:
    1. Read a txt file of fasta sequences  
    2. save names in one variable and sequences in another
    3. search the matches, mismatches, gaps and unknowns for each pairwise 
       alignment and sum the score
    4. Sum up the identities for each pairwise alignment
    5. Write an output file with the names of the individuals, the identity
       score and the alignment score for each pair
    
Input: 
    - text file with DNA sequences
    - [optional] text file with weight of matches, mismatches, gaps and unknowns 
      following this structure:
        match=[number]
        mismatch=[number]
        gap=[number]
        unknown=[number]

Output: one txt file with all posible pairwise alignments, their identity scores
        and their alignment score
    
Usage: 

    python MSA.py fasta_file weight_parameters[optional] output_file
    
"""

import os 
import sys
import pandas as pd


#####Basic control check#####

print("Checking input files: ...")

#for input length
if len(sys.argv) == 3: #when user inputs the output file
    fasta_user = sys.argv[1] #first argument is fasta file
    print("The fasta file that you provided is: ", sys.argv[1])
    output_file = sys.argv[2] #output file is the user's name    
    print("The output of this process is going to be stored the file named: ", sys.argv[2])
    match=5 
    mismatch=-4
    gap=-6
    unknown=-1 
elif len(sys.argv) == 2: #when user inputs just one file
    fasta_user=sys.argv[1] #first argument is fasta file
    print("The fasta file that you provided is: ", sys.argv[1])
    output_file = "output.txt" #output file is aminoacid_sequences.txt
    print("The output of this process is going to be stored the file named: output.txt")
    match=5 
    mismatch=-4
    gap=-6
    unknown=-1 
elif len(sys.argv) == 1: #when user just types the name of the code
    print("Please provide a fasta file and use this code like this: python MSA.py fasta_file weight_parameters[optional] output_file")
    sys.exit(1) #if not then exit
elif len(sys.argv) == 4: #when the user types more than expected
    fasta_user = sys.argv[1] #first argument is fasta file
    print("The fasta file that you provided is: ", sys.argv[1])
    weight_user = sys.argv[2] #second argument is the weight file
    print("The weights for the MSA are stored in the file named: ", sys.argv[2])
    if not os.path.exists(weight_user): #Does the file exist?
        print("The file you uploaded does not exist. Please provide an existing file.")
        sys.exit(1) #if not then exit
    elif not os.path.getsize(weight_user)>0: #Does the file contain values?
        print("The weight file you uploaded does not contain any value. Please provide a valid file following this structure: match=[number] \n mismatch=[number] \n gap=[number] \n unknown=[number]")
        sys.exit(1) #if not then exit
    with open (weight_user, "r") as weights: #open the file containing the weights
        for parameter in weights: #for each parameter in the file
            if parameter.startswith("match="): #if the line starts with match
                match = int(parameter[6:].strip()) #save the value in the variable match
            elif parameter.startswith("mismatch="): #if the line starts with mismatch
                mismatch = int(parameter[9:].strip()) #save the value in the variable match
            elif parameter.startswith("gap="): #if the line starts with gap
                gap = int(parameter[4:].strip()) #save the value in the variable match
            elif parameter.startswith("unknown="): #if the line starts with unknown
                unknown = int(parameter[8:].strip()) #save the value in the variable match
    output_file = sys.argv[3] #output file is the user's name    
    print("The output of this process is going to be stored the file named: ", sys.argv[3])
elif len(sys.argv) >= 5: #when the user types more than expected
    print("There are too many files, did you typed something extra? Please use this code like this: python MSA.py fasta_file weight_parameters[optional] output_file")
    sys.exit(1) #if not then exit

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

print("Checking finished")
#####Commands#####

#if theres a file then store the results in this, otherwuise use default as follows

print("Evaluating Multiple Sequence Alignment ...")

data = open(fasta_user,"r") #open the sequences file as data
sequences = [] #emty variable to store sequences
names=[] #empty variable to store names
for line in data: #for each line in the sequences file
    if not line.startswith(">"): #if is not a header
        sequences.append(line.strip()) #save the sequence
    else: #if its the header
        names.append(line[1:].strip()) #save the name
        

individual_scores=[] #empty matrix to save the individual scores  
for num1, seq1 in enumerate(sequences):  # iterate through the sequences
    s1 = names[num1] #save the name of the sequence with the index 
    for num2, seq2 in enumerate(sequences): #iterate through all the sequences
        s2 = names[num2] #save the name of the sequence with the index
        score = 0  # reset the score
        identity = 0 #reset the identity score
        for i in range(len(seq1)):  # for each position in the sequence until the end
            if seq1[i] in "ACGT" and seq2[i] in "ACGT":  # if both chars are bases
                if seq1[i] == seq2[i]:  # if there is a match add match score and add 1 to identity score
                    score += match
                    identity += 1
                else:  # if there is a mismatch add mismatch score
                    score += mismatch
            elif seq1[i] == "-" or seq2[i] == "-":  # if one of the chars is a gap add gap score
                if seq1[i] == seq2[i] == "-":  # if both chars are gaps add 0
                    score += 0
                else:  # if only one char is a gap add gap score
                    score += gap
            elif seq1[i] == "?" or seq2[i] == "?":  # if one of the chars is unknown add unknown score
                score += unknown
        identity = str(round(identity/len(seq1)*100, 1))+"%" #save identity as the percentage of identity
        inds = [str(s1) ,  str(s2) , str(identity), str(score)] #save a new variable 'inds' with both names, the identity scores and the overall score
        individual_scores.append(inds) #add 'inds' to the individual scores variable
#print(individual_scores)
#len(individual_scores[0])

print("Writing output file ...")
#make individual scores variable a datafrmae with column names as follows:
individual_scores_df=pd.DataFrame(individual_scores, columns=["SampleA","SampleB","IdentityScore","Score"])
#save the dataframe to an output tab separated without the column indexes
individual_scores_df.to_csv(output_file, sep="\t", index=False)

print("Process complete")

#close files
data.close()
