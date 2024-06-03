#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Name:FastaParser.py
Author:Maria Laura Mahecha Escobar
Description: Reads the genetic data #2 file split it into mtDNA and Y chromosome files
Institution:Lund University-Programming in python course
email:ma5668ma-s@student.lu.se
Date:2023y10m19d
Version:Python3.0 or more

User defined functions
- 

Procedure:
    1. Read a txt file sequences  
    2. search coincidences between names 
    3. save mtDNA sequence and y chromosome sequence for each name
    4. Write mtDNA sequences into one file
    5. Write Y chromosome sequences into another file
    
Input: text file with DNA sequences
Output: two fasta files containing mtDNA sequences (mtDNA.fna) and Y chromosome 
        sequences (Ychromosome.fna)
    
Usage: 

    python FastaParser.py text_file output_fasta_file_mtDNA[optional] output_fasta_file_Ychromosome[optional]
    
"""

#save the names of the individuals to match
individuals_names=["Princess Irene", "Prince Fred", "Nicolas II Romanov", "Alexandra Romanov", "Olga Romanov", "Tatiana Romanov", "Maria Romanov", "Alexei Romanov", "Suspected body of Anastasia Romanov", "Anastasia1", "Anastasia2", "Anastasia3", "Anastasia4", "Anastasia5", "Farmer's daughter", "Farmer’s grandson", "Grigori Rasputin"]

genetic_data = open("GeneticData - 2.txt", "r") #open the txt file
mt_dna = open("mtDNA.fna", "w") #open the mtDNA file to write
y_chromosome = open("Ychromosome.fna", "w") #open the Y chromosome file to write
for line in genetic_data: #for each line in the text file provided
    if any(name in line for name in individuals_names): #if the line contains a string that matches any of the names provided
        ind_name = line.strip() #save the name 
        if ind_name.startswith(">"): #if the line starts with >
            ind_name = ind_name[1:] #save the name without the >
    elif ind_name is not None: #if the name is not empty
        if line.startswith("mtDNA"): #if the line starts with mtDNA
            mt_seq = next(genetic_data).strip() #save all the lines that are sequence
            mt_dna.write(">" + ind_name + "\n" + mt_seq + "\n") #write in the mtDNA file a > the name and a new line with the sequence
        if line.startswith("Y chromosome"): #if the line starts with Y chromosome
            y_seq = next(genetic_data).strip() #save all lines that are sequence
            y_chromosome.write(">" + ind_name + "\n" + y_seq + "\n") #write in the Y chromosome file a > the name and a new line with the sequence

#close files
genetic_data.close()
mt_dna.close()
y_chromosome.close()
