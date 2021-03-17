#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Title:  Progenepredict: A python-based tool for de novo gene prediction in prokaryotic genomes

Date:   14.03.2021

Author: Lennart Hohmann

Description:
    Short: The program takes a fasta file containing one prokaryotic genome as the input file and predicts genes based
    on the presence of a Shine-Dargarno sequence before an open reading-frame (ORF).
    Long: The program creates the complement reverse strand to the leading strand that was input in fasta format.
    Then, it iterates over both strand, creating the three distinct reading-frame codons at each step.
    For each reading-frame the program searches for a start codon that is preceded by a Shine-Dalgarno sequence
    (assessed using a local alignment algorithm). If such a codon is found the program then searches for the next stop
    codon in that frame. After finding the stop codon, the length of the ORF is checked and if it fulfills the minimum
    length criterion, the predicted gene with associated information is written to the output files.

List of functions:
    1. get_revseq(sequence)
    2. SD_check(sequence, index, SD_seq)
    3. matrix(a, b, match_score=3, gap_cost=2)
    4. traceback(H, b, b_='', old_i=0)
    5. smith_waterman(a, b, match_score=3, gap_cost=2)
    6. detect_genes(input_file, min_length, start_codons, stop_codons)

List of "non standard" modules:
    1. NumPy (v1.19.4)

Procedure:
    1. Handle the user input, store it in the designated variables and print the settings the program is run with.
    2. The input file with the genome in fasta format is read and the leading strand saved.
    3. The reverse complement strand is created based on the leading strand.
    4. Go over the leading strand in steps of three and created the codons for the reading frames 1,2 and 3 during
       each iteration.
    5. For each reading frame codon it is checked whether it is 1. a start codon, is preceded by a Shine-Dalgarno
       sequence and there is currently no active start codon (currently only searching for an end codon) or 2. a
       stop codon and there is currently an active start codon.
    6. When a stop codon is found after a start codon preceded by a Shine-Dalgarno seq., it checks if the minimum ORF
       length criterion is fulfilled and if yes writes the predicted gene with associated information to the output file.
    7. Steps 4, 5, and 6 are then executed for the reverse strand as well.

Usage: progenepredict.py [-h] -g [-o] [-ml] [-sc  [...]] [-ec  [...]] [-sd] [-me]

Upcoming:
    1. User interface
    2. Simplify usage by automatically translating lowercase input of start/stop codons and Shine-Dalgarno sequences to
       uppercase
    3. Check that the input SD seq is not longer than 20 characters (longer than pre_sc seq)
"""

###########################################################

# import modules
import argparse
import itertools
import numpy as np

###########################################################

# defining functions

'''
Function(s): 
    get_revseq(sequence)
Goal: 
    Create the reverse complement DNA sequence to the input DNA sequence. 
    Used to create the reverse complement strand to the leading strand.
Key arguments:
    sequence(str) - A DNA sequence
'''

def get_revseq(sequence):
    reverse_seq = sequence[::-1].translate(str.maketrans("ATGC", "TACG"))
    return reverse_seq

'''
Function(s): 
    SD_check(sequence, index, SD_seq)
Goal: 
    Check whether a Shine-Dalgarno sequence is present in the 20 nucleotides before the start codon or not by applying 
    a local alignment algorithm and then checking if the number of errors exceeds the specified threshold to consider 
    a match valid.
Key arguments:
    sequence(str) - A DNA sequence
    index(int) - An index position in the previously specified DNA sequence
    SD_seq(str) - A shorter DNA sequence, in this application corresponding to the user specified Shine-Dalgarno sequence
'''

def SD_check(sequence, index, SD_seq):
    pre_sc = sequence[index-20:index] # only the 20 nucleotides before the start codon are checked
    sm_results = smith_waterman(pre_sc, SD_seq, match_score=3, gap_cost=2) # get the 1. start and 2. end alignment position coordiantes of the SD_seq in the pre_sc and 3. the string with how the SD_seq aligns to the pre_sc sequence (e.g. A-GGAGG)
    SD_alignment = sm_results[2] # this is how the SD sequence aligns to the pre_sc sequence (it may include gaps and/or deletions, e.g. A-GGAGG)
    match_in_pre_sc = pre_sc[sm_results[0]:sm_results[1]]  # store the sequence that the SD_alignment matches to in the pre_sc sequence (e.g. ATGGAGG)
    error_count = sum(1 for a, b in zip(match_in_pre_sc, SD_alignment) if a != b) + abs(len(match_in_pre_sc) - len(SD_alignment)) # count the number of positions at which the two strings differ from each other, if the strings differ in length that is considered as well
    if error_count > allowed_errors: # check against the setting of how many errors are allowed to occur for a match to be considered valid
        return False # not present
    else:
        return True # SD seq present

'''
Function(s): 
    matrix(a, b, match_score=3, gap_cost=2)
    traceback(H, b, b_='', old_i=0)
    smith_waterman(a, b, match_score=3, gap_cost=2)
Goal: 
    Perform a local alignment based on the Smith Waterman algorithm to find where the Shine-Dalgarno sequence (including gaps and deletions) best aligns to the 20 nucleotides before a start codon and retrieve the exact position in these nucleotides by using indeces of that alignment. The final function smith_waterman() returns the start and end coordinates of where the alignment starts and ends within these 20 nucleotides.
Key arguments:
    1. matrix(a, b, match_score=3, gap_cost=2)
    a(str) - the pre start codon sequence (longer than b)
    b(str) - the Shine-Dalgarno sequence (shorter than a)
    match_score=3(int) - the score that is added if a match is found and also corresponds to the penalty for mismatches (on that case it gets substracted)
    gap_cost=2(int) - the score that is subtracted for gaps
    
    2. traceback(H, b, b_='', old_i=0)
    H(numpy array) - numpy array of type int with the dimensions of the sequences a and b
    b(str) - the Shine-Dalgarno sequence
    b_=''(str) - string in which the best alignment of the Shine Dalgarno sequence is gradually saved
    old_i=0(int) - row index of the previously highest value in the matrix
    
    3. smith_waterman(a, b, match_score=3, gap_cost=2)
    a(str)- the pre start codon sequence (longer than b)
    b(str) - the Shine-Dalgarno sequence (shorter than a)
    match_score=3(int) - the score that is added if a match is found and also corresponds to the penalty for mismatches (on that case it gets substracted)
    gap_cost=2(int) - the score that is subtracted for gaps
    
References and resources:
    https://numpy.org/doc/stable/reference/
    https://tiefenauer.github.io/blog/smith-waterman/
'''

# Creating the scoring matrix
def matrix(a, b, match_score=3, gap_cost=2):
    H = np.zeros((len(a) + 1, len(b) + 1), int) # Returns a new array of type int filled with zeros
    for i, j in itertools.product(range(1, H.shape[0]), range(1, H.shape[1])): # go though the array. i = row; j = column
        match = H[i - 1, j - 1] + (match_score if a[i - 1] == b[j - 1] else - match_score) # diagonal step (-1 for i and j because it takes the value of the cell 1 row and 1 column before), if the two sequences match at that position the match score is added. If not the missmatch is penalized by substracting the match score value
        delete = H[i - 1, j] - gap_cost # step to the side (previsous row position i-1, same column position), so the gap cost is substracted
        insert = H[i, j - 1] - gap_cost # step down (previous column position j-1, same row position), so the gap cost is substracted
        H[i, j] = max(match, delete, insert, 0) # check which of these options is the best and set the current array position to that value
    return H

# Backtracing the alignment in the matrix
def traceback(H, b, b_='', old_i=0):
    # it is necessary to do the flip step because the indices of the LAST occurrence of H.max() has to be returned by np.argmax()
    H_flip = np.flip(np.flip(H, 0), 1) # 1st flip turns the array on its head; 2nd flip then turns it like a page in a book so the bottom right corner is the top left corner now
    # H_flip = np.flip(H) # could also flip like this, but the previous line is clearer of whats going on
    i_, j_ = np.unravel_index(H_flip.argmax(), H_flip.shape) # np.unravel_index(indices, shape) converts a flat index into a tuple of coordinate arrays. So np.unravel_index(index of the maximum value, array dimensions) returns the coordinates of the maximum value in the matrix
    i, j = np.subtract(H.shape, (i_ + 1, j_ + 1))  # i,j is the index pair of the last occurrence of the H.max() value
    if H[i, j] == 0: # the function is continuously returned until the traceback finds a 0
        #print(b_,j)
        return b_, j # returns the final b_ which is the way b aligns to a (including gaps and deletions)
    b_ = b[j - 1] + '-' + b_ if old_i - i > 1 else b[j - 1] + b_ # in the backtracing the letters are added one by one to b_, but if the new i (row of the current max value) is more than 1 step away from the row of the previsous maximum value, it means that the next highest value was not in the diagonal field and a gap is introduced
    #print(b_)
    return traceback(H[0:i, 0:j], b, b_, i) # the considered matrix gets smaller and smaller, always keeping the maximum value in the bottom right corner (i,j are the indeces of that value), the alignment is gradually saved in b_

# Calculating the start and end index
def smith_waterman(a, b, match_score=3, gap_cost=2):
    a, b = a.upper(), b.upper() # not necessarily needed here as a,b will be upper already but useful when utilizing the function in other contexts
    H = matrix(a, b, match_score, gap_cost) # create the matrix
    b_, pos = traceback(H, b) # backtrace
    return pos, pos + len(b_), b_ # for this application i need to return the start, end, AND the b_ alignment for later comparison (to evaluate the match of the SD seq to the pre start codon seq)

'''
Function(s): 
    detect_genes(input_file, min_length, start_codons, stop_codons)
Goal: 
    Read in the input file with the genome in fasta format and save the leading strand. Based on the leading strand the
    reverse complement strand is created. Then the function loops over the leading strand in steps of three and created 
    the codons for reading frames 1,2 and 3 each step. For each reading frame codon it is checked whether it is 1. a 
    start codon, is preceded by a Shine-Dalgarno sequence and there is currently no active start codon (currently 
    only searching for an end codon) or 2. a stop codon and there is currently an active startcodon. When a stop codon 
    is found it checks if the minimum ORF length criterion is fulfilled and if yes writes the predicted gene with 
    associated information to the output file. The same process then happends for the reverse strand.
    
Key arguments:
    input_file(str) - An input file containing one prokaryotic genome in fasta format
    min_length(int) - A minimum ORF length that a gene has to have to be predicted as one
    start_codons(set) - A set with the start codons
    stop_codons(set) - A set with the stop codons
'''

# finding the ORFs
def detect_genes(input_file, min_length, start_codons, stop_codons, fasta_output):
    # first the input file is processed and the forward strand saved
    with open(input_file,'r') as genome:
        forward_strand = ''
        for line in genome:
            if line.startswith('>'): # there should be only one line with the name of the genome
                genome_name = line.strip().lstrip('>') # not further used for now
            else:
                forward_strand += line.strip() # if the sequence is split over multiple lines it will be concatenated
        # based on the forward strand the reverse complement strand is created
        reverse_strand = get_revseq(forward_strand)
    # variables for the three frames
    f1_startflag = False
    f2_startflag = False
    f3_startflag = False
    geneid_count = 0
    with open(output_file, "w") as output:
        output.write('#{}\t{}\t{}\t{}\t{}\t{}\n'.format('gene_id','start','end','readingframe','strand','start_type')) # the column names
        # open a second output file if the user wants the predicted egene sequences stored in fasta format
        if fasta_output:
            fastafile = open("predictedgenes.fasta","w")
        else:
            pass
        # create the three reading frames for the forward strand
        for i in range(0,len(forward_strand),3):  # create the index value, go in steps of three. i will be 0,3,6,etc.
            # for frame one
            f1_codon = forward_strand[i:i+3] # create the current codon
            if len(f1_codon) == 3: # check that is not the end of the sequence
                if f1_codon in start_codons and SD_check(forward_strand, i, SD_seq) and not f1_startflag: # check if it is a start codon preceeded by a SD sequence and there is no current start codon (position of codon matters, can be start codon or simply coding for the given amino acid)
                    f1_startflag = True # a start codon was found
                    f1_sc = f1_codon # save the start codon
                    f1_startindex = i # save the index where the ORF starts
                elif f1_startflag and f1_codon in stop_codons: # check if the codon is a stop codon
                    f1_endindex = i+3 # save the index where the ORF ends (i is the start of the stopcodon, so i+3)
                    f1_startflag = False # when the stop codon is found the startflag is set to false
                    if (f1_endindex - f1_startindex) >= min_length: # check if it qualifies as an ORF with the minimum length
                        geneid_count += 1 # variable to name the found genes
                        output.write('gene_{}\t{}\t{}\t{}\t{}\t{}\n'.format(geneid_count, f1_startindex, f1_endindex, 1, '+',f1_sc)) # + here because in forward strand, 1 because readingframe 1
                        if fasta_output:
                            fastafile.write('>gene_{}\tstart={}\tend={}\treadingframe={}\tstrand={}\tstart_type={}\n{}\n'.format(geneid_count, f1_startindex, f1_endindex, 1, '+',f1_sc,forward_strand[f1_startindex:f1_endindex]))
                        else:
                            pass
                    else:
                        pass
                else:
                    pass
            else:
                pass
            # for frame two
            f2_codon = forward_strand[i+1:i+4]  # create the current codon
            if len(f2_codon) == 3:  # check that is not the end of the sequence
                if f2_codon in start_codons and SD_check(forward_strand,i,SD_seq) and not f2_startflag:  # check if it is a start codon preceeded by a SD sequence and there is no current start codon
                    f2_startflag = True  # a start codon was found
                    f2_sc = f2_codon  # save the start codon
                    f2_startindex = i+1  # save the index where the ORF starts
                elif f2_startflag and f2_codon in stop_codons:  # check if the codon is a stop codon
                    f2_endindex = i+4  # save the index where the ORF ends (i is the start of the stopcodon, so i+3)
                    f2_startflag = False  # when the stop codon is found the startflag is set to false
                    if (f2_endindex - f2_startindex) >= min_length:  # check if it qualifies as an ORF with the minimum length
                        geneid_count += 1  # variable to name the found genes
                        output.write('gene_{}\t{}\t{}\t{}\t{}\t{}\n'.format(geneid_count,f2_startindex,f2_endindex,2,'+',f2_sc))  # + here because in forward strand, 2 because readingframe 2
                        if fasta_output:
                            fastafile.write('>gene_{}\tstart={}\tend={}\treadingframe={}\tstrand={}\tstart_type={}\n{}\n'.format(geneid_count,f2_startindex,f2_endindex,2,'+',f2_sc,forward_strand[f2_startindex:f2_endindex]))
                        else:
                            pass
                    else:
                        pass
                else:
                    pass
            else:
                pass
            # for frame 3
            f3_codon = forward_strand[i+2:i+5]  # create the current codon
            if len(f3_codon) == 3:  # check that is not the end of the sequence
                if f3_codon in start_codons and SD_check(forward_strand,i,SD_seq) and not f3_startflag:  # check if it is a start codon preceeded by a SD sequence and there is no current start codon
                    f3_startflag = True  # a start codon was found
                    f3_sc = f3_codon  # save the start codon
                    f3_startindex = i+2  # save the index where the ORF starts
                elif f3_startflag and f3_codon in stop_codons:  # check if the codon is a stop codon
                    f3_endindex = i+5  # save the index where the ORF ends (i is the start of the stopcodon, so i+3)
                    f3_startflag = False  # when the stop codon is found the startflag is set to false
                    if (f3_endindex - f3_startindex) >= min_length:  # check if it qualifies as an ORF with the minimum length
                        geneid_count += 1  # variable to name the found genes
                        output.write('gene_{}\t{}\t{}\t{}\t{}\t{}\n'.format(geneid_count, f3_startindex, f3_endindex, 3,'+',f3_sc))  # + here because in forward strand, 3 because readingframe 3
                        if fasta_output:
                            fastafile.write('>gene_{}\tstart={}\tend={}\treadingframe={}\tstrand={}\tstart_type={}\n{}\n'.format(geneid_count, f3_startindex, f3_endindex, 3,'+',f3_sc,forward_strand[f3_startindex:f3_endindex]))
                        else:
                            pass
                    else:
                        pass
                else:
                    pass
            else:
                pass

        # create the three reading frames for the reverse strand
        for i in range(0, len(reverse_strand),3):  # create the index value, go in steps of three. i will be 0,3,6,etc.
            # for frame one
            f1_codon = reverse_strand[i:i+3]  # create the current codon
            if len(f1_codon) == 3:  # check that is not the end of the sequence
                if f1_codon in start_codons and SD_check(reverse_strand, i,SD_seq) and not f1_startflag:  # check if it is a start codon preceeded by a SD sequence and there is no current start codon
                    f1_startflag = True  # a start codon was found
                    f1_sc = f1_codon  # save the start codon
                    f1_startindex = i  # save the index where the ORF starts
                elif f1_startflag and f1_codon in stop_codons:  # check if the codon is a stop codon
                    f1_endindex = i+3  # save the index where the ORF ends (i is the start of the stopcodon, so i+3)
                    f1_startflag = False  # when the stop codon is found the startflag is set to false
                    if (f1_endindex - f1_startindex) >= min_length:  # check if it qualifies as an ORF with the minimum length
                        geneid_count += 1  # variable to name the found genes
                        output.write('gene_{}\t{}\t{}\t{}\t{}\t{}\n'.format(geneid_count, f1_startindex, f1_endindex, 1,'-',f1_sc))  # - here because in forward strand, 1 because readingframe 1
                        if fasta_output:
                            fastafile.write('>gene_{}\tstart={}\tend={}\treadingframe={}\tstrand={}\tstart_type={}\n{}\n'.format(geneid_count, f1_startindex, f1_endindex, 1,'-',f1_sc,reverse_strand[f1_startindex:f1_endindex]))
                        else:
                            pass
                    else:
                        pass
                else:
                    pass
            else:
                pass
            # for frame two
            f2_codon = reverse_strand[i+1:i+4]  # create the current codon
            if len(f2_codon) == 3:  # check that is not the end of the sequence
                if f2_codon in start_codons and SD_check(reverse_strand, i,SD_seq) and not f2_startflag:  # check if it is a start codon preceeded by a SD sequence and there is no current start codon (open ORF)
                    # make a flag so that start codons should not be considred until the next stop codon is found
                    f2_startflag = True  # a start codon was found
                    f2_sc = f2_codon  # save the start codon
                    f2_startindex = i+1  # save the index where the ORF starts
                elif f2_startflag and f2_codon in stop_codons:  # check if the codon is a stop codon
                    f2_endindex = i+4  # save the index where the ORF ends (i is the start of the stopcodon, so i+3)
                    f2_startflag = False  # when the stop codon is found the startflag is set to false
                    if (f2_endindex - f2_startindex) >= min_length:  # check if it qualifies as an ORF with the minimum length
                        geneid_count += 1  # variable to name the found genes
                        output.write('gene_{}\t{}\t{}\t{}\t{}\t{}\n'.format(geneid_count, f2_startindex, f2_endindex, 2,'-',f2_sc))  # - here because in forward strand, 2 because readingframe 2
                        if fasta_output:
                            fastafile.write('>gene_{}\tstart={}\tend={}\treadingframe={}\tstrand={}\tstart_type={}\n{}\n'.format(geneid_count, f2_startindex, f2_endindex, 2,'-',f2_sc,reverse_strand[f2_startindex:f2_endindex]))
                        else:
                            pass
                    else:
                        pass
                else:
                    pass
            else:
                pass
            # for frame 3
            f3_codon = reverse_strand[i+2:i+5]  # create the current codon
            if len(f3_codon) == 3:  # check that is not the end of the sequence
                if f3_codon in start_codons and SD_check(reverse_strand, i,SD_seq) and not f3_startflag:  # check if it is a start codon preceeded by a SD sequence and there is no current start codon (open ORF)
                    # make a flag so that start codons should not be considred until the next stop codon is found
                    f3_startflag = True  # a start codon was found
                    f3_sc = f3_codon  # save the start codon
                    f3_startindex = i+2  # save the index where the ORF starts
                elif f3_startflag and f3_codon in stop_codons:  # check if the codon is a stop codon
                    f3_endindex = i+5  # save the index where the ORF ends (i is the start of the stopcodon, so i+3)
                    f3_startflag = False  # when the stop codon is found the startflag is set to false
                    if (f3_endindex - f3_startindex) >= min_length:  # check if it qualifies as an ORF with the minimum length
                        geneid_count += 1  # variable to name the found genes
                        output.write('gene_{}\t{}\t{}\t{}\t{}\t{}\n'.format(geneid_count, f3_startindex, f3_endindex, 3,'-',f3_sc))  # - here because in forward strand, 3 because readingframe 3
                        if fasta_output:
                            fastafile.write('>gene_{}\tstart={}\tend={}\treadingframe={}\tstrand={}\tstart_type={}\n{}\n'.format(geneid_count, f3_startindex, f3_endindex, 3,'-',f3_sc,reverse_strand[f3_startindex:f3_endindex]))
                        else:
                            pass
                    else:
                        pass
                else:
                    pass
            else:
                pass
        # if fasta_ouput == True then the file was opened and has to be closed
        if fasta_output:
            fastafile.close()
        else:
            pass


###########################################################

# running code

# handling user input
parser = argparse.ArgumentParser(description="This program takes a fasta file containing one prokaryotic genome and predicts genes based on the presence of a Shine-Dargarno sequence before an open reading frame (ORF).")
parser.add_argument('-g', '--genomefile', help="The file containing the genome sequence in fasta format.", metavar='', required=True)
parser.add_argument('-o', '--outputfile', help="The desired name of the output file containing the predicted genes. (default=progenepredict.txt)", default='progenepredict.txt', metavar='')
parser.add_argument('-ml', '--minORFlength', help="The minimum ORF length of predicted genes. (default=300)", metavar='', default=300, type = int)
parser.add_argument('-sc', '--startcodons', help="The possible start codons. Input like: -sc ATG GTG TTG ; (default=ATG)", nargs='+', metavar='', default=["ATG"], type=str)
parser.add_argument('-ec', '--stopcodons', help="The possible stop codons. Input like: -ec TAA TAG TGA ; (default=TAA TAG TGA)", nargs='+', metavar='', default=["TAA","TAG","TGA"], type=str)
parser.add_argument('-sd', '--shinedalgarnosequence', help="The Shine-Dalgarno sequence. (default=AGGAGG)", metavar='', default="AGGAGG", type=str)
parser.add_argument('-me', '--maxalignerrors', help="The maximal number of errors that are allowed in the pairwise alignemnt when determining whether a Shine-Dalgarno sequence is present or not. (default=1)", metavar='', default=1, type=int)
parser.add_argument('-f', '--fastaoutput', help="Output a file with the DNA sequences of the predicted genes in fasta format", action='store_true')

args = parser.parse_args()

# saving user input in specific variables
start_codons = set(args.startcodons)    # convert the list with start codons into a set
stop_codons = set(args.stopcodons)      # convert the list with stop codons into a set
input_file = args.genomefile
output_file = args.outputfile
min_length = args.minORFlength
SD_seq = args.shinedalgarnosequence
allowed_errors = args.maxalignerrors
fasta_output = args.fastaoutput

# print out the settings the program is run with
print('The program is run with the settings:\ninput genome file:\t{}\noutput file name:\t{}\nminimum ORF length:\t{}\nShine-Dalgarno-seq.:\t{}\nstart codons:\t\t{}\nstop codons:\t\t{}\nalign. error threshold:\t{}\nouput fasta:\t\t{}\n'.format(input_file, output_file, min_length, SD_seq, start_codons, stop_codons, allowed_errors, fasta_output))

if __name__ == '__main__': # only execute the code when it is run as a program
    detect_genes(input_file, min_length, start_codons, stop_codons, fasta_output)
