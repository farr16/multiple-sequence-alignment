#!/usr/bin/python

################################################
# Lab7 CS423 starter file
# Put your name here
# Fall 2015
################################################

import random, math

######################################################################
# Determine the score of the optimal global alignment of two sequences
# copy this function from lab 5
# Just calculates the score (not the alignment; do not need to
# print out tables; just return the optimal score)
######################################################################
def globalAlignmentScore(s1, s2):

    # Scoring system
    MATCH = 5
    MISMATCH = -4
    GAP = -6

    # Set table size
    NUM_ROWS = len(s2) + 1
    NUM_COLS = len(s1) + 1

    # Create table and fill it with zeroes
    costs = createTable(NUM_ROWS, NUM_COLS, 0)

    # Fill the top row of the costs table with values decrementing by 6
    val = 0
    for i in range(0, NUM_COLS):
        costs[0][i] = val
        val -= 6

    # Fill the left column of the costs table with values decrementing
    # by 6
    val = -6
    for j in range(1, NUM_ROWS):
        costs[j][0] = val
        val -= 6

    # Fill in each of the rows with alignment costs, starting at the
    # top row going left to right
    for y in range(1, NUM_ROWS):
        for x in range(1, NUM_COLS):
            # Calculate costs of a gap for the top and bottom sequences
            valTop = costs[y-1][x] + GAP
            valLeft = costs[y][x-1] + GAP

            valDag = 0
            # Increase score for a match, or decrease for a mismatch
            if s1[x-1] == s2[y-1]:
                valDag = costs[y-1][x-1] + MATCH
            else:
                valDag = costs[y-1][x-1] + MISMATCH

            # Calculate the maximum value and set that one as our cost
            val = max(valTop, valLeft, valDag)
            costs[y][x] = val

    # Return optimal score (lower right-hand cell in the table)
    return costs[NUM_ROWS - 1][NUM_COLS - 1]


#################################################################
# Generate a random DNA sequence with given length and nucleotide
# probabilities
# complete this function
# A previous lab exercise may be helpful
#################################################################
def generateRandomSequence(length, A_prob, C_prob, G_prob, T_prob):
    return ""	# delete this stub


##################################################################
# Genenerate a random DNA sequence comparable in nucleotide
# distribution as the given input sequence
# complete this function
##################################################################
def generateComparableRandomSequence(s):
    # calculate nucleotide frequencies and call generateRandomSequence function
    # to create the random DNA sequence similar in nucleotide composition and same length
    return ""	# delete this stub


####################################################################
# For the two input sequences s1 and s2, calculate the distance score D where
# D is definted as 100.0*(-ln(S_norm)). 
#
# S_norm = (S_global - S_rand) / (S_iden - S_rand)
#
# S_global is the optimal global alignment score between s1 and s2.
#
# S_iden is the average of the global alignment scores of s1 aligned
# with s1 and s2 aligned with s2.
#
# S_rand is the average of 1000 global alignment scores between
# sequences similar in composition as s1 and s2.
# complete this function
###################################################################
def calculateDistanceScore(s1, s2):
    return 0	# delete this stub

####################################################################
# Create a 2D table with the given number of rows and columns
# and fills all entries with value given as a parameter
# (function completed for you)
####################################################################
def createTable(numRows, numCols, value):
     table = []
     row = 0
     # create 2D table initialized with value
     while (row < numRows):
          table.append([])    
          col = 0
          while (col < numCols):
               table[row].append(value)
               col = col + 1
          row = row + 1
     return table


########################################################
### End of functions ###################################
########################################################

seq1 = "GCGTGCATGCATGCCGGG"
seq2 = "ATTAGACGGATATTCGAT"
score = globalAlignmentScore(seq1, seq2)
print(score)

### DNA sequences
##seq1 = "CGATAGTGCTATATCTAGCGCCGTCTAGATGCATTATACGATATCG"
##seq2 = "AACGACATGGCTCGTGCTATTACGCGCGAATATCC"
##seq3 = "ATAGTGCTATACTCGTGCTATTCTAGATGCCGCGATATAT"
##seq4 = "GGATAGGCTATATCTAGCGCGTCTAGATGCATTTACGATATC"
##seq5 = "TACGACATGCGCTCGTGCATATTAGCGCGCGATATATCG"
##
### calculate global alignment scores for all ten pairs
##print("Global alignment scores")
##print("seq1 and seq2: " + str(globalAlignmentScore(seq1, seq2)))
### finish for other nine pairs
##
### calculate distance score for all ten pairs
##print("")
##print("Distance scores")
##print("seq1 and seq2: " + str(calculateDistanceScore(seq1, seq2)))
### finish for other nine pairs


