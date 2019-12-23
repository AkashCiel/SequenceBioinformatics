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
referenceGenome: TextIO = open("PindelReferenceGenome.txt")
mainFile: TextIO = open("SamFileSmall.txt")

# Read reference genome
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
readLength = 50
insertSize = 100
maxDeletionSize = 200

for i in range(totalReadPairs):
    allMinMaxSubstrings.append([])
    # If left read is mapped and right is not, start looking from 3' end of left read
    if np.char.find(refGenome[0], readPairs[i][0]) > -1 and np.char.find(refGenome[0], readPairs[i][1]) == -1:
        anchorPoint = np.char.find(refGenome[0], readPairs[i][0]) + readLength
        minMax3prime = minmaxSubstrings(refGenome[0][anchorPoint:anchorPoint + 2*insertSize], readPairs[i][1], reverse=False)
        if minMax3prime != []:
            threePrimeAnchor = np.char.find(refGenome[0][anchorPoint:anchorPoint + 2*insertSize], minMax3prime[0])
            minMax5prime = minmaxSubstrings(refGenome[0][anchorPoint + threePrimeAnchor:anchorPoint + threePrimeAnchor + readLength + maxDeletionSize],
                                            readPairs[i][1], reverse=True)
            allMinMaxSubstrings[i].append(1)
            allMinMaxSubstrings[i].append(minMax3prime)
            allMinMaxSubstrings[i].append(minMax5prime)
        else:
            allMinMaxSubstrings[i].append(1)
            allMinMaxSubstrings[i].append([])
            allMinMaxSubstrings[i].append([])

    # If right read is mapped and left if not, look upto the 3' end of right read
    elif np.char.find(refGenome[0], readPairs[i][0]) == -1 and np.char.find(refGenome[0], readPairs[i][1]) > -1:
        anchorPoint = np.char.find(refGenome[0], readPairs[i][1])
        minMax3prime = minmaxSubstrings(refGenome[0][max(0, anchorPoint - 2*insertSize):anchorPoint], readPairs[i][0], reverse=True)
        if minMax3prime != []:
            threePrimeAnchor = np.char.find(refGenome[0][max(0, anchorPoint - 2*insertSize):anchorPoint][::-1], minMax3prime[0][::-1])
            minMax5prime = minmaxSubstrings(refGenome[0][max(0, anchorPoint - threePrimeAnchor - readLength - maxDeletionSize)
                                                         :anchorPoint - threePrimeAnchor], readPairs[i][0], reverse=True)
            allMinMaxSubstrings[i].append(0)
            allMinMaxSubstrings[i].append(minMax3prime)
            allMinMaxSubstrings[i].append(minMax5prime)
        else:
            allMinMaxSubstrings[i].append(0)
            allMinMaxSubstrings[i].append([])
            allMinMaxSubstrings[i].append([])
    # If both reads are mapped or none is, there is nothing to be done
    elif np.char.find(refGenome[0], readPairs[i][1]) > -1 and np.char.find(refGenome[0], readPairs[i][0]) > -1:
        allMinMaxSubstrings[i] == None
    elif np.char.find(refGenome[0], readPairs[i][1]) == -1 and np.char.find(refGenome[0], readPairs[i][0]) == -1:
        allMinMaxSubstrings[i] == None

breakPoints = []
for i in range(totalReadPairs):
    if allMinMaxSubstrings[i] != []:
        if allMinMaxSubstrings[i][1] != [] and allMinMaxSubstrings[i][2] != []:
            allPossibleReconstructions = []
            pairIndex = allMinMaxSubstrings[i][0]
            for x in range(len(allMinMaxSubstrings[i][1])):
                for y in range(len(allMinMaxSubstrings[i][2])):
                    reconstruction = allMinMaxSubstrings[i][1][x] + allMinMaxSubstrings[i][2][y]
                    if reconstruction == readPairs[i][pairIndex]:
                        breakPoint01 = np.char.find(refGenome[0], allMinMaxSubstrings[i][1][x]) + len(allMinMaxSubstrings[i][1][x])
                        breakPoint02 = np.char.find(refGenome[0], allMinMaxSubstrings[i][2][y]) + 0
                        breakPoints.append([i, pairIndex, breakPoint01, breakPoint02])
                    allPossibleReconstructions.append(reconstruction)
            if readPairs[i][pairIndex] in allPossibleReconstructions:
                print(i, readPairs[i][pairIndex])
                print(allMinMaxSubstrings[i])
for item in breakPoints:
    print(item)

