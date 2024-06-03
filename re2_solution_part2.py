#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Title: Toolbox RE2 - Part 2 solution
Created on: 2022-10-14
Author: Arthur Boffelli Castro

Description:
    This script contains a possible solution for the Running Exercise II -
        Part 2.
    This solution is only an example of how to achieve the expected result, the
        exercise may be solved in many different ways.

    The script contains three independent programs:
        dna2aa      - Translates DNA sequences from a fasta or fastq file to
                      amino acid sequences.
        amino_count - Counts the total number of each amino acids found in a
                      fasta or fastq file.
        barcode     - removes three possible 8 nucleotide barcodes from
                      sequences in a fastq or fasta file.

    All three programs require a fasta or fastq file to run, the output files
        are optional, if not provided the program generates an output file with
        a default name.

Procedure:
    1. Parse all arguments from the command line using argparse.
    2. Check if the input files exist and are in the corresponded formats.
    3. Run the block of code according to the program chosen by the user.
        a. dna2aa - Checks if the sequence lines contain only DNA bases. If not
            the user is questioned if they would like to continue. The sequence
            lines are then read in triplets (codons) and translated to their
            respective amino acid or stop codon, that are store in a dictionary.
            Codons containing non DNA bases are translated to "X".
        b. amino_count - Reads only the sequence lines and using a dictionary to
            store the frequency of each amino acid, outputs a list of amino
            acids and their respective frequencies throughout the whole file.
        c. barcode - Isolate the first 8 characters from the sequence and
            check if the 8-mer is one of the three possible barcodes. Next,
            check the last 8 characters in the same procedure. If a barcode is
            found the barcode is removed and the sequence is added to the output
            file. If none of the barcodes are found, the read is added to the
            undetermined output.

Usage:
    python re2_solution_part2.py [-h] (--fasta IN_FASTA | --fastq IN_FASTQ) [-o OUTFILE] [-u UNDETERMINED] PROGRAM

Examples:
    python re2_solution_part2.py dna2aa --fasta fasta_example.fasta -o amino_out.fasta
    python re2_solution_part2.py --fastq barcode.fastq barcode
"""

import argparse
import sys
import re


#%% Argparse section.
usage = '''Tool box containing three different programs (all three programs \
require a fasta or fastq input file to run):\n'''
programs = """\
\tdna2aa\t\tTranslation of DNA sequences to amino acid sequences.
\tamino_count\tCount the number of each amino acid found in a amino acid fasta \
file.
\tbarcode\t\tRemove barcodes from the beginning or end of DNA sequences in a \
fastq file. 
"""
# Create the argparse functionality, the description will be displayed in the
# help option. formatter_class is used to read the string in raw format,
# allowing the use of \n and \t.
parser = argparse.ArgumentParser(description=usage + programs,
                                 formatter_class=argparse.RawTextHelpFormatter)

# Create a positional argument, this argument does not have a flag.
parser.add_argument('PROGRAM', type=str, action="store",
                    help="One of the three possible programs")

# Create a group of arguments that are required, this allows us to have two
# arguments and one of them is required.
group = parser.add_mutually_exclusive_group(required=True)

# Add argument with flag --fasta to the required group.
group.add_argument("--fasta", dest='infile', type=str, metavar="IN_FASTA",
                   help="Input file in fasta format")

# Add argument with flag --fastq to the required group.
group.add_argument('--fastq', dest='in_fastq', type=str, metavar="IN_FASTQ",
                   help="Input file in fastq format")

# Add optional argument with flag to output file.
parser.add_argument('-o', '--outfile', dest='outfile', type=str, default='',
                    help="Output file name [default <PROGRAM>_output.txt]")

# Add optional argument with flag to undetermined barcode sequences.
parser.add_argument('-u', '--undetermined', dest='undetermined', type=str,
                    default="undetermined.txt",
                    help="Additional output file for the barcode program which "
                         "includes sequences without barcodes [default "
                         "undetermined.txt]")

# Parse the arguments and save it in a variable, this creates a class object.
args = parser.parse_args()

#%% General functions


def fasta_qc(fasta_file):
    """
    Function to check if the input file is in fasta format.

    Parameters
    ----------
    fasta_file : str
        Path for single of multi line fasta file.

    Raises
    ------
    TypeError
        If format is not correct.

    Returns
    -------
    bool
        True if the format is correct.
        None if the format is not correct.

    """
    try:
        with open(fasta_file, 'r') as fasta:
            # Check if the first line starts with >.
            first_line = fasta.readline()
            if first_line.startswith('>'):
                return True
            else:
                raise TypeError
    # Add error message for file not found
    except FileNotFoundError as not_found:
        print("The file {} was not found!".format(not_found.filename))

    # Add error message for file not in fasta format
    except TypeError as type_error:
        print("Input file must be in fasta format!")


def fastq_qc(fastq_file):
    """
    Function to check if the input file is in fastq format.

    Parameters
    ----------
    fastq_file : string
        Path for fastq file.

    Raises
    ------
    TypeError
        If format is not correct.

    Returns
    -------
    bool
        True if the format is correct.
        None if the format is not correct.

    """
    try:
        with open(fastq_file, 'r') as fastq:
            # Check if the first line starts with @.
            first_line = fastq.readline()
            if first_line.startswith('@'):
                return True
            else:
                raise TypeError
    # Add error message if file not found
    except FileNotFoundError as not_found:
        print("The file {} was not found!".format(not_found.filename))

    # Add error message for file not in fastq format
    except TypeError as type_error:
        print("Input file must be in fastq format!")


def fasta_bases_qc(input_file):
    """
    Function to check if the sequence lines in a fasta file contains only 
    DNA bases.
    If different letters are found in the sequence lines, the user is asked to
    continue or not (with "y" or "n"). The default answer is "n", which will
    stop the program.

    Parameters
    ----------
    input_file : str
        Path for fasta file.

    Returns
    -------
    bool
        True if only DNA bases are found or if the user answers "y" to continue.
        Exit the program if answer is "n".
    """
    with open(input_file) as fasta:
        # Start an empty set
        bases_set = set()
        for fasta_line in fasta:
            # We are interested only in the sequence lines, so ignore the
            # headers
            if not fasta_line.startswith(">"):
                # Transform the line to lower case and create a set, the set
                # will keep only unique letters.
                line_set = set(fasta_line.strip().lower())
                # Join the set created to the previous set, so we check all the
                # unique characters in all sequence lines.
                bases_set = bases_set | line_set

        # Join the set into a string.
        fasta_bases = ''.join(sorted(list(bases_set)))

        # Using regular expression, we check if our string contains any
        # character that is not a DNA nucleotide. This function returns True if
        # any other character except [a,t,c,g] is found.
        if re.search(r"[^acgt]", fasta_bases):

            # Ask the user if they want to continue even though there are non
            # DNA characters in the sequence.
            proceed = input("The file contain non-DNA bases. Would you "
                            "like to continue (codons containing non-"
                            "DNA letters will be translated into X? "
                            "y/[n] ").lower()
            if proceed == 'y':  # User chooses to continue.
                return True

            else:  # User cancels the run.
                print("Cancelled!")
                sys.exit(0)

        else:  # No other character was found, only DNA characters.
            return True


def fastq_bases_qc(input_file):
    """
    Function to check if the sequence lines in a fastq file contains only
    DNA bases.
    If different letters are found in the sequence lines, the user is asked to
    continue or not (with "y" or "n"). The default answer is "n", which will
    stop the program.

    Parameters
    ----------
    input_file : str
        Path for fastq file.

    Returns
    -------
    bool
        True if only DNA bases are found or if the user answers "y" to continue.
        Exit the program if answer is "n".
    """
    with open(input_file) as fastq:
        # Start an empty set
        bases_set = set()
        for fastq_line in fastq:
            # Keep track of the four lines of the fastq file.
            if fastq_line.startswith("@"):
                sequence = next(fastq)
                spacer = next(fastq)
                quality = next(fastq)
                # Transform the line to lower case and create a set, the set
                # will keep only unique letters.
                line_set = set(sequence.strip().lower())
                # Join the set created to the previous set, so we check all the
                # unique characters in all sequence lines.
                bases_set = bases_set | line_set

        # Join the set into a string.
        fastq_bases = ''.join(sorted(list(bases_set)))

        # Using regular expression, we check if our string contains any
        # character that is not a DNA nucleotide. This function returns True if
        # any other character except [a,t,c,g] is found.
        if re.search(r"[^acgt]", fastq_bases):

            # Ask the user if they want to continue even though there are non
            # DNA characters in the sequence.
            proceed = input("The file contain non-DNA bases. Would you "
                            "like to continue (codons containing non-"
                            "DNA letters will be translated into X? "
                            "y/[n] ")

            if proceed == 'y':  # User chooses to continue.
                return True
            else:  # User cancels the run.
                print("Cancelled!")
                sys.exit(0)
        else:  # No other character was found, only DNA characters.
            return True


def translation(dna_sequence):
    """
    Function to translate a DNA sequence into an amino acid sequence.

    Parameters
    ----------
    dna_sequence: str
        Sequence containing DNA nucleotides.

    Returns
    -------
    aa_sequence: str
        Sequence of amino acids translated from the DNA sequence.
    """
    # Make the whole sequence uppercase and transform Ts to Us.
    dna_sequence = dna_sequence.upper().replace("T", "U")
    # Initialize the amino acid string.
    aa_sequence = ''

    # When looping through the DNA sequence, we subtract two from the total
    # length to make sure we are getting only triples. If the last bases do not
    # form a codon, they are ignored. We loop in steps of three.
    for i in range(0, len(dna_sequence) - 2, 3):
        # Index the initial letter and the following two.
        codon = dna_sequence[i:i + 3]
        # Retrieve the amino acid for the respective codon from the codon
        # dictionary. Here we set a default, if the key does not exist, we use
        # "X" instead.
        aa_sequence += codon_table.get(codon, "X")
    return aa_sequence


#%% Preparation of variables

# Part one - dna2aa ############################################################
# Create a codon dictionary using the 4 nucleotides and all amino acids in the
# right order as strings.
bases = "UCAG"
amino_acids = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'

# First we create a list with all possible codons, using list comprehension.
codons = [a + b + c for a in bases for b in bases for c in bases]

# Zip the codon list with the amino acids string and transform it into a
# dictionary.
codon_table = dict(zip(codons, amino_acids))

# Part two - amino_count #######################################################
# Create a set with all amino acids using the previous dictionary values.
amino_acids_set = set(codon_table.values())

# Create a dictionary to count the amino acids, set all values to zero.
count_dict = dict(zip(sorted(codon_table.values()),
                      [0 for i in range(len(codon_table))]))

# Remove the stop codon key from the dictionary, since it is not a amino acid.
count_dict.pop("*")

# Add the key X to count the letters that are not amino acids.
count_dict["X"] = 0

# Part three - barcode #########################################################
barcodes = ["TATCCTCT", "GTAAGGAG", "TCTCTCCG"]

# Input and output files #######################################################
# Retrieve the positional argument that has which program will be used.
prog_type = args.PROGRAM.lower()

# Check if the program name is one of the possible three, if not print a help
# message and exit.
if prog_type not in ["dna2aa", "amino_count", "barcode"]:
    print("Command not recognised. Please use one of the following options:")
    print(programs)
    sys.exit(1)

# Assign all arguments to their respective variables.
infasta = args.infile
infastq = args.in_fastq
outfile = args.outfile
undetermined = args.undetermined

# If output file is empty, join the name of the program chosen to a default
# output name.
if not outfile:
    outfile = prog_type + "_output.txt"

#%% Running code for fasta
# Part of the code that deals with fasta files.

# If the infasta variable has a value, and it passes on the fasta_qc (if the
# file exists, and it is a fasta file).
if infasta and fasta_qc(infasta):

    # Part one - dna2aa.
    if prog_type == "dna2aa":

        # Continue only if the bases_qc return True (either by having only DNA
        # bases or the user continuing after the message).
        if fasta_bases_qc(infasta):

            # Open the input file and the output file.
            with open(infasta, 'r') as input_fasta, open(outfile, 
                                                         'w') as out_file:

                # To be able to deal with multi line fasta file we create a
                # while loop. We initiate an empty sequence and read the lines
                # one by one.
                sequence = ''
                while True:
                    line = input_fasta.readline()

                    # When we have a header or the file ends (empty line), check
                    # if there is a sequence.
                    if line.startswith(">") or not line:
                        if sequence:  # Skip if it is the first header.
                            # Translate the DNA sequence and print to the output
                            # file.
                            amino_seq = translation(sequence)
                            print(amino_seq, file=out_file)
                            # Break the while loop if reached the end of the
                            # file
                            if not line:
                                break
                        # Empty the sequence variable to start a new one.
                        sequence = ''
                        # Save the header and print to the output file.
                        header = line.strip()
                        print(header, file=out_file)
                    else:  # sequence line.
                        # Sum to the existing sequence.
                        sequence += line.strip()

                print("Done!")

    # Part two - amino_count.
    elif prog_type == "amino_count":

        # Open the input file and the output file.
        with open(infasta, 'r') as input_fasta, open(outfile, 'w') as out_file:
            for line in input_fasta:
                # In this case we are only interested on the sequence lines.
                if not line.startswith(">"):
                    # Loop through each character of the amino acid sequence.
                    for aa in line.strip().upper():
                        # If the letter is an amino acid, sum up in the count.
                        if aa in amino_acids_set:
                            count_dict[aa] += 1
                        else:  # Non-amino acid letter, sum up the letter X.
                            count_dict["X"] += 1

            # Print the count dictionary to the output file.
            for aa in count_dict:
                print(f"{aa}\t{count_dict[aa]}", file=out_file)

            print("Done!")

    # Part three - barcode.
    # barcode is meant to be used in a fastq file, however the program
    # accepts it in fasta format as well.
    elif prog_type == "barcode":

        # Open the input file and the two output files.
        with open(infasta, 'r') as input_fasta, open(
                outfile, 'w') as out_file, open(undetermined, 'w') as undet:

            # To be able to deal with multi line fasta file we create a
            # while loop. We initiate an empty sequence and read the lines
            # one by one.
            sequence = ''
            while True:
                line = input_fasta.readline()

                # When we have a header or the file ends (empty line), check
                # if there is a sequence.
                if line.startswith(">") or not line:
                    if sequence:  # Skip if it is the first header.
                        # Copy the sequence to a new variable.
                        trimmed = sequence

                        # The barcodes are only present in the start or end
                        # of the sequences. Check the first 8 characters and
                        # the last 8 characters of the sequence.
                        # Note that I am using all if statements instead of
                        # elif, this ensures that both start and end of
                        # the same sequence will be checked (removing a
                        # possible double barcode)
                        if sequence[:8] in barcodes:
                            trimmed = trimmed[8:]
                        if sequence[-8:] in barcodes:
                            trimmed = trimmed[:-8]

                        # Check if the trimmed sequence has a different
                        # length and save to the trimmed file or
                        # undetermined file if nothing changed.
                        if len(sequence) != len(trimmed):
                            print(f"{header}\n{trimmed}", file=out_file)
                        else:  # length is the same.
                            print(
                                f"{header}\n{trimmed}", file=undet)

                    # Break the while loop if reached the end of the
                    # file
                    if not line:
                        break

                    # Empty the sequence variable to start a new one and
                    # save the header.
                    sequence = ''
                    header = line.strip()

                else:  # sequence line
                    # Sum to the existing sequence.
                    sequence += line.strip().upper()
            print("Done!")

#%% Running code for fastq

# If the infastq variable has a value, and it passes on the fastq_qc (if the
# file exists, and it is a fastq file).
elif infastq and fastq_qc(infastq):

    # Part one - dna2aa.
    # dna2aa is meant to be used in a fasta file, however the program
    # accepts it in fastq format as well, and outputs a fasta file.
    if prog_type == "dna2aa":

        # Continue only if the bases_qc return True (either by having only DNA
        # bases or the user continuing after the message).
        if fastq_bases_qc(infastq):
            # Open the input file and the output file.
            with open(infastq, 'r') as input_fastq, open(
                    outfile, 'w') as out_file:

                # Since fastq files are the worst type of files for
                # bioinformaticians, the best way to keep track of the lines is
                # storing the four lines from the beginning.
                for line in input_fastq:
                    if line.startswith("@"):
                        # From the first line, store all four lines using the
                        # next function.
                        header = line.strip()
                        sequence = next(input_fastq).strip()
                        spacer = next(input_fastq).strip()
                        quality = next(input_fastq).strip()

                        # Translate the sequence line.
                        amino_seq = translation(sequence)

                        # Print the header and the amino acid sequence in fasta
                        # format.
                        print('>' + header[1:], file=out_file)
                        print(amino_seq, file=out_file)

                print("Done!")

    # Part two - amino_count.
    # amino_count is meant to be used in a fasta file, however the program
    # accepts it in fastq format as well.
    elif prog_type == "amino_count":

        # Open the input file and the output file.
        with open(infastq, 'r') as input_fastq, open(outfile, 'w') as out_file:

            # Keep track of the four lines in the fastq file since the
            # beginning.
            for line in input_fastq:
                if line.startswith("@"):
                    sequence = next(input_fastq)
                    spacer = next(input_fastq)
                    quality = next(input_fastq)

                    # Loop through all characters in the sequence line.
                    for aa in sequence.strip().upper():
                        # If the letter is an amino acid, sum up in the count.
                        if aa in amino_acids_set:
                            count_dict[aa] += 1
                        else:  # Non-amino acid letter, sum up the letter X.
                            count_dict["X"] += 1

            # Print the count dictionary to the output file.
            for aa in count_dict:
                print(f"{aa}\t{count_dict[aa]}", file=out_file)

            print("Done!")

    # Part three - barcode.
    elif prog_type == "barcode":

        # Open the input file and the two output files.
        with open(infastq, 'r') as input_fastq, open(
                outfile, 'w') as out_file, open(undetermined, 'w') as undet:

            # Keep track of the four lines in the fastq file since the
            # beginning.
            for line in input_fastq:
                if line.startswith("@"):
                    header = line.strip()
                    sequence = next(input_fastq).strip().upper()
                    spacer = next(input_fastq).strip()
                    quality = next(input_fastq).strip()
                    # Copy the sequence to a new variable.
                    trimmed = sequence

                    # The barcodes are only present in the start or end
                    # of the sequences. Check the first 8 characters and
                    # the last 8 characters of the sequence.
                    # Note that I am using all if statements instead of
                    # elif, this ensures that both start and end of
                    # the same sequence will be checked (removing a
                    # possible double barcode)
                    if sequence[:8] in barcodes:
                        trimmed = trimmed[8:]
                        quality = quality[8:]
                    if sequence[-8:] in barcodes:
                        trimmed = trimmed[:-8]
                        quality = quality[:-8]

                    # Check if the trimmed sequence has a different
                    # length and save the four lines to the trimmed file or
                    # undetermined file if nothing changed.
                    if len(sequence) != len(trimmed):
                        print(f"{header}\n{trimmed}\n{spacer}\n{quality}",
                              file=out_file)
                    else:
                        print(f"{header}\n{trimmed}\n{spacer}\n{quality}",
                              file=undet)

            print("Done!")

else:  # input file failed the qc_function.
    # Quit the program.
    sys.exit(1)
    # NOTE: The use of 0 or 1 in the exit is not python specific. A non-zero
    # exit means that there was some kind of error/issue/problem during the
    # execution of the code. A zero exit means a successful exit.

