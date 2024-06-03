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
import os

packages=["sys", "pandas"] #Make a list of the modules needed for this script to run

for i in packages: #Install the modules if they are not installed in the computer of the person running the code
    try:
        __import__(i)
    except:
        os.system("pip install "+i)

import os
import sys
import pandas as pd


print("Running script to add the BLAST hit description to the fasta sequence header: ", sys.argv[0])

#File input from the user
fasta_user=sys.argv[1]
print("The fasta file that you provided is: ", sys.argv[1])
blast_user=sys.argv[2]
print("The file containing the hitDescription and sequence ID is: ", sys.argv[2])
output_file=sys.argv[3]
print("The output of this process is going to be in the file named: ", sys.argv[3],"\n")

print("Checking input files: ...")



###Basic control check###
#for fasta file

if not os.path.exists(fasta_user): #Does the file exist?
    print("The fasta file you uploaded does not exist. Please provide an existing file.")
    sys.exit(1) #if not then exit
elif not os.path.getsize(fasta_user)>0: #Does the file contain values?
    print("The fasta file you uploades does not contain any sequence. Please provide a sequence")
    sys.exit(1) #if not then exit

with open (fasta_user, "r") as fasta_file:
    first_line_fasta=fasta_file.readline() #Does the first line of the fasta file starts with something different of ">"
    if not first_line_fasta.startswith (">"):
        raise TypeError("It seems you uploaded a non supported type of file, maybe a Fastq? Please check and provide a fasta format file")
        sys.exit(1) #if yes then exit
    else:
        print("Fasta file seems to be correct")


#For BLAST file
if not os.path.exists(blast_user): #Does the file exist?
    print("The BLAST file you uploaded does not exist. Please provide an existing file.")
    sys.exit(1) #if not then exit
elif not os.path.getsize(blast_user)>0:  #Does the file contain values?
    print("The BLAST file you uploaded does not contain any sequence. Please provide a sequence")
    sys.exit(1) #if not then exit
else:
    print("BLAST file seems to be correct")

print("Checking finished.","\n","Running...")



###Comands###


#Open blast file to turn into a dictionary
blast_file=open(blast_user, "r")
#Store queryname and hit description as a dictionary only when hit description is != null
#uses pandas
#read the blast file as a csv delimited by tabs and counting first row as column names
blast_csv=pd.read_csv(blast_file, sep='\t', lineterminator='\n', header=(0))
#make the csv a dataframe to work with
blast_df=pd.DataFrame(blast_csv)
#extract ONLY the query name and the hit descriptions 
blast_df_cols=blast_df[['#queryName', 'hitDescription']]
#extract the rows that have all values (not include rows with null values)
blast_df_cols=blast_df_cols.dropna()
#Make a dictionary that sets the queryname as key () and the hit description as values []
blast_dic=blast_df_cols.set_index("#queryName")["hitDescription"].to_dict()
#blast_dic["2_g"]
#len(blast_dic)

if len(blast_dic)!=0: #Check how many sequences have are left after extracting the 'null' values
    print(sys.argv[2], "has", len(blast_dic), "sequences with hitDescription different from 'null'")
else:
    print(sys.argv[2], "has no sequences with hitDescription different from 'null'")
    exit()


#open a new fasta to store results
new_fasta_file = open (output_file, "w")

print("Adding hitDescription to fasta sequences")

#make a variable to store temporary results
new_line=None
#open the original fasta file
with open (fasta_user, "r") as fasta_file:
    for line in fasta_file: #for each line in the fasta file
        if line.startswith (">"):
            header = line[1:] #save the header withput the > chr
            header = header.split() #save it with the tabs
            if header[0] in blast_dic: #if the first position (ID) appears in the blast dictionary
                new_line = line.strip() + "\t" + "protein="+"".join(blast_dic[header[0]]) #store in the new line the line from the fasta file (without removing any chr, that's strip() for) and the value from the dictionary
                new_fasta_file.write (new_line + "\n") #writes this line in the output file
        elif new_line is not None: #To add the next line in the output file, check if there was something in the new_line variable
            new_fasta_file.write(line) #print the line of the iteration (it does not start with ">")
            new_line=None #set new_line variable to zero

print("Process complete")


#closing all files
fasta_file.close()
blast_file.close()
new_fasta_file.close()