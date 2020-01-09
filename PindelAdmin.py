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
totalReadPairs = 1635

# Read DNA sequences from text
# Import text files in Python
referenceGenome: TextIO = open("PindelReferenceLarge.txt")
mainFile: TextIO = open("SamFileLarge.txt")

# Read reference genome
refGenomeSeq = referenceGenome.readlines()
refGenomeSTR = ''
for i in range(0, len(refGenomeSeq)):
    refGenomeSeq[i] = refGenomeSeq[i].rstrip('\n')
    refGenomeSTR += refGenomeSeq[i]

refGenome = np.array([refGenomeSTR])
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
readLength = 100
insertSize = 200
maxDeletionSize = 400

# Looking for deletions in the reference genome
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
            #print(np.char.find(refGenome[0], minMax3prime[0]) + len(minMax3prime[1]), readPairs[i][1], allMinMaxSubstrings[i])
        else:
            allMinMaxSubstrings[i].append(1)
            allMinMaxSubstrings[i].append([])
            allMinMaxSubstrings[i].append([])

    # If right read is mapped and left is not, look upto the 3' end of right read
    elif np.char.find(refGenome[0], readPairs[i][0]) == -1 and np.char.find(refGenome[0], readPairs[i][1]) > -1:
        anchorPoint = np.char.find(refGenome[0], readPairs[i][1])
        minMax5prime = minmaxSubstrings(refGenome[0][max(0, anchorPoint - 2*insertSize):anchorPoint], readPairs[i][0], reverse=True)
        if minMax5prime != []:
            fivePrimeAnchor = np.char.find(refGenome[0][max(0, anchorPoint - 2*insertSize):anchorPoint][::-1], minMax5prime[0][::-1])
            minMax3prime = minmaxSubstrings(refGenome[0][max(0, anchorPoint - fivePrimeAnchor - readLength - maxDeletionSize)
                                                         :anchorPoint - fivePrimeAnchor], readPairs[i][0], reverse=False)
            allMinMaxSubstrings[i].append(0)
            allMinMaxSubstrings[i].append(minMax3prime)
            allMinMaxSubstrings[i].append(minMax5prime)
            #if minMax3prime != []:
            #    print(np.char.find(refGenome[0], minMax3prime[0]) + len(minMax3prime[1]), readPairs[i][0], allMinMaxSubstrings[i])
        else:
            allMinMaxSubstrings[i].append(0)
            allMinMaxSubstrings[i].append([])
            allMinMaxSubstrings[i].append([])
    # If both reads are mapped or none is, there is nothing to be done
    elif np.char.find(refGenome[0], readPairs[i][1]) > -1 and np.char.find(refGenome[0], readPairs[i][0]) > -1:
        allMinMaxSubstrings[i] == None
    elif np.char.find(refGenome[0], readPairs[i][1]) == -1 and np.char.find(refGenome[0], readPairs[i][0]) == -1:
        allMinMaxSubstrings[i] == None

breakPointsDeletions = []
for i in range(totalReadPairs):
    if allMinMaxSubstrings[i] != []:
        if allMinMaxSubstrings[i][1] != [] and allMinMaxSubstrings[i][2] != []:
            if len(allMinMaxSubstrings[i][1][1]) + len(allMinMaxSubstrings[i][2][1]) >= readLength:
                pairIndex = allMinMaxSubstrings[i][0]
                leftBreakPoint = np.char.find(refGenome[0], allMinMaxSubstrings[i][1][1]) + len(allMinMaxSubstrings[i][1][1])
                rightBreakPoint = np.char.find(refGenome[0], allMinMaxSubstrings[i][2][1]) + 0
                breakPointsDeletions.append([leftBreakPoint, rightBreakPoint])

# Looking for insertion events in the reference genome
allMinMaxSubstrings = []
for i in range(totalReadPairs):
    allMinMaxSubstrings.append([])
    # If left read is mapped and right is not, start looking from 3' end of left read
    if np.char.find(refGenome[0], readPairs[i][0]) > -1 and np.char.find(refGenome[0], readPairs[i][1]) == -1:
        anchorPoint = np.char.find(refGenome[0], readPairs[i][0]) + readLength
        minMax3prime = minmaxSubstrings(refGenome[0][anchorPoint:anchorPoint + 2*insertSize], readPairs[i][1], reverse=False)
        if minMax3prime != []:
            threePrimeAnchor = np.char.find(refGenome[0][anchorPoint:anchorPoint + 2*insertSize], minMax3prime[0])
            minMax5prime = minmaxSubstrings(refGenome[0][anchorPoint + threePrimeAnchor:anchorPoint + threePrimeAnchor + readLength - 1],
                                            readPairs[i][1], reverse=True)
            allMinMaxSubstrings[i].append(1)
            allMinMaxSubstrings[i].append(minMax3prime)
            allMinMaxSubstrings[i].append(minMax5prime)
            #print(np.char.find(refGenome[0], minMax3prime[0]) + len(minMax3prime[1]), readPairs[i][1], allMinMaxSubstrings[i])
        else:
            allMinMaxSubstrings[i].append(1)
            allMinMaxSubstrings[i].append([])
            allMinMaxSubstrings[i].append([])

    # If right read is mapped and left is not, look upto the 3' end of right read
    elif np.char.find(refGenome[0], readPairs[i][0]) == -1 and np.char.find(refGenome[0], readPairs[i][1]) > -1:
        anchorPoint = np.char.find(refGenome[0], readPairs[i][1])
        minMax5prime = minmaxSubstrings(refGenome[0][max(0, anchorPoint - 2*insertSize):anchorPoint], readPairs[i][0], reverse=True)
        if minMax5prime != []:
            fivePrimeAnchor = np.char.find(refGenome[0][max(0, anchorPoint - 2*insertSize):anchorPoint][::-1], minMax5prime[0][::-1])
            minMax3prime = minmaxSubstrings(refGenome[0][max(0, anchorPoint - fivePrimeAnchor - readLength + 1)
                                                         :anchorPoint - fivePrimeAnchor], readPairs[i][0], reverse=False)
            allMinMaxSubstrings[i].append(0)
            allMinMaxSubstrings[i].append(minMax3prime)
            allMinMaxSubstrings[i].append(minMax5prime)
            #if minMax3prime != []:
            #    print(np.char.find(refGenome[0], minMax3prime[0]) + len(minMax3prime[1]), readPairs[i][0], allMinMaxSubstrings[i])
        else:
            allMinMaxSubstrings[i].append(0)
            allMinMaxSubstrings[i].append([])
            allMinMaxSubstrings[i].append([])
    # If both reads are mapped or none is, there is nothing to be done
    elif np.char.find(refGenome[0], readPairs[i][1]) > -1 and np.char.find(refGenome[0], readPairs[i][0]) > -1:
        allMinMaxSubstrings[i] == None
    elif np.char.find(refGenome[0], readPairs[i][1]) == -1 and np.char.find(refGenome[0], readPairs[i][0]) == -1:
        allMinMaxSubstrings[i] == None

breakPointsInsertions = []
for i in range(totalReadPairs):
    if allMinMaxSubstrings[i] != []:
        if allMinMaxSubstrings[i][1] != [] and allMinMaxSubstrings[i][2] != []:
            breakPoint01 = np.char.find(refGenome[0], allMinMaxSubstrings[i][1][1]) + len(allMinMaxSubstrings[i][1][1])
            insertionLength = readLength - (len(allMinMaxSubstrings[i][1][1]) + len(allMinMaxSubstrings[i][2][1]))
            breakPointsInsertions.append([breakPoint01, insertionLength])

for item in breakPointsInsertions:
    if breakPointsInsertions.count(item) < 2:
        breakPointsInsertions.remove(item)

print("Reporting deletions: ")
breakPointsDeletions.sort(key=itemgetter(0))
for item in breakPointsDeletions:
    print(item)

print("Reporting insertions: ")
breakPointsInsertions.sort(key=itemgetter(0))
for item in breakPointsInsertions:
    print(item)

