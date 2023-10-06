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
import pandas as pd


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




# =============================================================================
#LETS BACK TO THIS LATER
# #For each line in fasta_file search for lines begining with ">" and extract all
# #information until another the next ">" appears and store the result in vasiable
# #named id_seq 
# 
# seq_fas=[] #open a variable to store temporarly all the fasta sequences
# new_fasta=[] #open a new variable to store all the modified sequences in fasta format
# 
# for line in fasta_file:
#     if not line.startswith(">"):
#         seq_fas.append(line)
#         if seq_fas:
#             new_fasta.append(seq_fas); #Keep the previous sequence in new_fasta
#         seq_fas='' #erases the content in seq_fas
#     else:
#         seq_fas.append(line)
# print(new_fasta)
# 
# 
# #An alternative shorter way to reading FASTA files
# r = open('Fasta_example_3_seq.txt', 'r')
# seq = [] #Holds the final combined sequences
# temp=''
# 
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
# 
# print("Printing the concatenated file:")
# print(seq)
# =============================================================================

#MARIA LAURA TENEMOS U DICCIONARIO PARA TRABAJAR DESPUÉS!!
#Mamasita, el trabajo que se te viene limpiando este codigo es monumental. Buena suerte belleza
#Y por si alguien se pregunta, si, hablo conmigomisma en tercera persona cuando escribo códigos. Me da moral cuando me siento frustrada, como ahora :)
#Y por si alguien llega a ver esto en el futuro, si, estoy siendo lenta e ineficiente pero que le hacemos, asi funciona mi mente :D bai por ahora
#Ahora si tenemos diccionario, no habia puesto el index pa'sortearlo, jejeje
#=============================================================================
#Store queryname and hit description as a dictionary only when hit description is != null
#uses pandas
# #read the blast file as a csv delimited by tabs and counting first row as column names
# blast_csv=pd.read_csv(blast_file, sep='\t', lineterminator='\n', header=(0))
# #make the csv a dataframe to work with
# blast_df=pd.DataFrame(blast_csv)
# #extract ONLY the query name and the hit descriptions 
# blast_df_cols=blast_df[['#queryName', 'hitDescription']]
# #extract the rows that have all values (not include rows with null values)
# blast_df_cols=blast_df_cols.dropna()
# #Seting the index from which values must be sorted in the dictionary
# blast_df_cols=blast_df_cols.set_index("#queryName")
# #Change column name 
# blast_df_cols.rename(columns = {'hitDescription':'protein'}, inplace = True)
# #make the new dataframe a dictionary
# blast_dic=blast_df_cols.to_dict("index")
# blast_dic["2_g"]
# =============================================================================







for lines in textfile.xreadlines():
    for eachkey in dict.keys():
        if eachkey in lines:
            print lines + " : " + dict[eachkey]
        else:
            continue








#Closing all files
fasta_file.close()
blast_file.close()
new_fasta_file.close()