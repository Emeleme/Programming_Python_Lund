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

print("Running script to add the BLAST hit description to the fasta sequence header: ", sys.argv[0])

fasta_file=sys.argv[1]
blast_file=sys.argv[2]
output_file=sys.argv[3]
    
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






#MARIA LAURA TENEMOS U DICCIONARIO PARA TRABAJAR DESPUÉS!!
#Mamasita, el trabajo que se te viene limpiando este codigo es monumental. Buena suerte belleza
#Y por si alguien se pregunta, si, hablo conmigomisma en tercera persona cuando escribo códigos. Me da moral cuando me siento frustrada, como ahora :)
#Y por si alguien llega a ver esto en el futuro, si, estoy siendo lenta e ineficiente pero que le hacemos, asi funciona mi mente :D bai por ahora
# =============================================================================

#Open blast file to turn into a dictionary
blast_file=open("five_seq.blastx.tab", "r")
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
blast_dic["2_g"]

#open a new fasta to store results
new_fasta_file = open ("output.txt", "w")

#THERE'S A PROBLEM WITH THIS CODE, I JOIN THE DICTIONARY BUT JUST ADDS THE WORD PROTEIN, I HAVE TO DELETE IT SOMEHOW
#Vamo'a intentar quitar los nombres de las columnas en el dataframe si mujer? solo a ver que pasa. Termina de comentar, le subimos a Github y seguimos 
#make a variable to store temporary results
new_line=None
#open the original fasta file
with open ("five_seq.fna", "r") as fasta_file:
    for line in fasta_file: #for each line in the fasta file
        if line.startswith (">"):
            header = line[1:] #save the header withput the > chr
            header = header.split() #save it with the tabs
            if header[0] in blast_dic: #if the first position (ID) appears in the blast dictionary
                new_line = line.strip() + "\t" + "".join(blast_dic[header[0]]) #store in the new line the line from the fasta file (without removing any chr, that's strip() for) and the value from the dictionary
                new_fasta_file.write (new_line + "\n") #writes this line in the output file
        elif new_line is not None: #To add the next line in the output file, check if there was something in the new_line variable
            new_fasta_file.write(line) #print the line of the iteration (it does not start with ">")
            new_line=None #set new_line variable to zero

#Closing all files
fasta_file.close()
blast_file.close()
new_fasta_file.close()