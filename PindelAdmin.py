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

from PindelFunctions import *

# Specify total number of sequences
totalReadPairs = 285

# Read DNA sequences from text
# Import text files in Python
reads01: TextIO = open("PindelReads01.txt")
reads02: TextIO = open("PindelReads02.txt")
referenceGenome: TextIO = open("PindelReferenceGenome.txt")
mainFile: TextIO = open("SamFileSmall.txt")

# Read files into text strings
Sequences = reads01.readlines()
leftEndReads = np.array([])
for i in range(0, len(Sequences)):
    Sequences[i] = Sequences[i].rstrip('\n')
    if i%2 == 1:
        leftEndReads = np.append(leftEndReads, Sequences[i])

Sequences = reads02.readlines()
rightEndReads = np.array([])
for i in range(0, len(Sequences)):
    Sequences[i] = Sequences[i].rstrip('\n')
    if i%2 == 1:
        rightEndReads = np.append(rightEndReads, Sequences[i])

refGenomeSeq = referenceGenome.readlines()
refGenome = np.array([])
for i in range(0, len(refGenomeSeq)):
    refGenomeSeq[i] = refGenomeSeq[i].rstrip('\n')
    refGenome = np.append(refGenome, refGenomeSeq[i])

# Read SAM
samReads = mainFile.readlines()
for i in range(0, len(samReads)):
    samReads[i] = samReads[i].rstrip('\n')
    samReads[i] = samReads[i].split("\t")

# Obtain read pairs from SAM
readPairs = []
for i in range(totalReadPairs):
    readPairs.append([])
    readPairs[i].append(samReads[2*i][9])
    readPairs[i].append(samReads[2*i + 1][9])

allMinMaxSubstrings = []

for i in range(totalReadPairs):
    allMinMaxSubstrings.append([])
    if np.char.find(refGenome[0], readPairs[i][0]) > -1 and np.char.find(refGenome[0], readPairs[i][1]) == -1:
        minMax5prime = minmaxSubstrings(refGenome[0], readPairs[i][1], reverse=False)
        minMax3prime = minmaxSubstrings(refGenome[0], readPairs[i][1], reverse=True)
        allMinMaxSubstrings[i].append(1)
        allMinMaxSubstrings[i].append(minMax5prime)
        allMinMaxSubstrings[i].append(minMax3prime)
    elif np.char.find(refGenome[0], readPairs[i][1]) > -1 and np.char.find(refGenome[0], readPairs[i][0]) == -1:
        minMax5prime = minmaxSubstrings(refGenome[0], readPairs[i][0], reverse=False)
        minMax3prime = minmaxSubstrings(refGenome[0], readPairs[i][0], reverse=True)
        allMinMaxSubstrings[i].append(0)
        allMinMaxSubstrings[i].append(minMax5prime)
        allMinMaxSubstrings[i].append(minMax3prime)
    elif np.char.find(refGenome[0], readPairs[i][1]) > -1 and np.char.find(refGenome[0], readPairs[i][0]) > -1:
        allMinMaxSubstrings[i] == None
    elif np.char.find(refGenome[0], readPairs[i][1]) == -1 and np.char.find(refGenome[0], readPairs[i][0]) == -1:
        allMinMaxSubstrings[i] == None



for i in range(totalReadPairs):
    if allMinMaxSubstrings[i] != []:
        if allMinMaxSubstrings[i][1] != [] and allMinMaxSubstrings[i][2] != []:
            allPossibleReconstructions = []
            for x in range(len(allMinMaxSubstrings[i][1])):
                for y in range(len(allMinMaxSubstrings[i][2])):
                    allPossibleReconstructions.append(allMinMaxSubstrings[i][1][x] + allMinMaxSubstrings[i][2][x])
            pairIndex = allMinMaxSubstrings[i][0]
            if readPairs[i][pairIndex] in allPossibleReconstructions:
                print(i, readPairs[i][pairIndex])
                print(allMinMaxSubstrings[i])

