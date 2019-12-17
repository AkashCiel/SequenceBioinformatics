from typing import TextIO
from operator import itemgetter
from random import seed
from random import randint
from collections import defaultdict
from copy import deepcopy
from random import randint
from decimal import *
import math
import numpy as np

# Convert number to base 2 string
def intoBase2String(number):
    binaryConversion = "{0:b}".format(number)
    return binaryConversion.zfill(12)

# Return reverse compliment of the string
def reverseCompliment(sequence):
    sequence = sequence[::-1]
    reverseSequence = ''
    for nucleotide in sequence:
        if nucleotide == "A": reverseSequence+="T"
        elif nucleotide == "T": reverseSequence+="A"
        elif nucleotide == "G": reverseSequence+="C"
        elif nucleotide == "C": reverseSequence+="G"
    return reverseSequence

# Find minimum and maximum unique substrings
def minmaxSubstrings(genome, pattern, reverse):
    if reverse == True:
        genome = genome[::-1]
        pattern = pattern[::-1]

    subStrings = []
    found = False
    while found == False:
        for i in range(1, len(pattern) + 1):
            subPattern = pattern[:i]
            if np.char.count(genome, subPattern) == 1:
                subStrings.append(subPattern)
            elif np.char.count(genome, subPattern) == 0:
                found = True
                break

    finalSubStrings = []
    if len(subStrings) > 0:
        finalSubStrings.append(subStrings[0][::((-1*(reverse == True)) + (1*(reverse == False)))])
        finalSubStrings.append(subStrings[len(subStrings) - 1][::((-1*(reverse == True)) + (1*(reverse == False)))])
    return finalSubStrings

# Detect a deletion event
def detectDeletion(genome, read01, read02, insertionSize, maxDeletionSize):
    # Locate read01 location on genome
    read01Anchor = np.char.find(genome, read01) + len(read01)
    # Find minimum maximum substrings from 3-prime end of second read
    minmax3prime = minmaxSubstrings(genome[read01Anchor:(read01Anchor + 2*insertionSize)], read02, reverse=False)
    if minmax3prime != []:
        # Anchor the beginning of substring defined above in the genome
        read02Anchor3prime = read01Anchor + np.char.find(genome[read01Anchor:(read01Anchor + 2*insertionSize)], minmax3prime[0])
        # Find minimum maximum substrings from 5-prime end of second read
        minmax5prime = minmaxSubstrings(genome[read02Anchor3prime:(read02Anchor3prime+maxDeletionSize+len(read02))], read02, reverse=True)
        # Reverese the substrings from 5-prime end, since read02 runs from 3-prime to 5-prime
        for i in range(len(minmax5prime)): minmax5prime[i] = minmax5prime[i][::-1]
        return [minmax3prime, minmax5prime]
    else:
        return []

