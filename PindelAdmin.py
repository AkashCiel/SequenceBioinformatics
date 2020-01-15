"""
This section contains code for the Pindel Administrator, the code block which orchestrates the entire process of reading
in the inputs, processing them in appropriate arrays, and then looking for insertions and deletion events. Kindly make
sure to input values for the appropriate variables as indicated below.
"""

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

#####################################################################
################ ENTER TOTAL NUMBER OF READ-PAIRS BELOW #############
#####################################################################

totalReadPairs = 1635

#################################################################################
################ PROVIDE NAME OF REFERENCE GENOME FILE (.TXT) BELOW #############
#################################################################################

referenceGenome: TextIO = open("PindelReferenceLarge.txt")

#################################################################################
#################### PROVIDE NAME OF SAM FILE (.TXT) BELOW ######################
#################################################################################

mainFile: TextIO = open("SamFileLarge.txt")

################################################################################################################
#################### PROVIDE READ LENGTH, INSERTION SIZE, AND MAXIMUM DELETION SIZE BELOW ######################
################################################################################################################

readLength = 100
insertSize = 200
maxDeletionSize = 400

####################################################################
#### Processing the input in appropriate arrays ####################
####################################################################

# Read reference genome

refGenomeSeq = referenceGenome.readlines()
refGenomeSTR = ''
for i in range(0, len(refGenomeSeq)):
    refGenomeSeq[i] = refGenomeSeq[i].rstrip('\n')
    refGenomeSTR += refGenomeSeq[i]

refGenome = np.array([refGenomeSTR])

# Read SAM File

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

"""
After processing the input sequences, the code
below attempts to map the read pairs sequentially against the reference genome. No action is needed or taken if both or
neither of the sequences from a read pair are mapped completely. If, however, only one of the sequences from a read-pair
can be mapped, it indicates the possibility of the other read coming from the site of an insertion/deletion event in the
reference genome. Using the function minmaxSubstrings (see documentation in PindelFunctions for details), we first look
for minimum and maximum unique substrings of the unmapped read from the end closer to the mapped read (the left end of
right read if the left read was mapped perfectly and right end of the left read if the right one was mapped), in a
smaller region of the reference genome. This new-found substring position (both the minimum and maximum unique substrings
will be positioned at the same spot in the reference genome), is used as an anchor to determine the region of reference
genome in which we again look for minimum-maximum unique substrings, from the opposite direction this time. The smaller
region of the reference genome we look for the substrings in is determined by the following factors:
    1. Are we looking for an insertion event, or deletion?
    2. Are we searching for min-max unique substrings from the left read, or the right one?
The exact formulae implemented in the code are as per the guidelines provided in the paper. If, using these tuples, we can
perfectly recreate the original unmapped read, it indicates a deletion event. If not, it indicates an insertion event.
The break-points are stored away for reference in either case.
"""

###########################################################
##### Looking for deletions in the reference genome #######
###########################################################

for i in range(totalReadPairs):
    allMinMaxSubstrings.append([])

    # If left read is mapped and right is not, start looking from the right end of left read (in the left-right direction)

    if np.char.find(refGenome[0], readPairs[i][0]) > -1 and np.char.find(refGenome[0], readPairs[i][1]) == -1:
        anchorPoint = np.char.find(refGenome[0], readPairs[i][0]) + readLength
        minMax3prime = minmaxSubstrings(refGenome[0][anchorPoint:anchorPoint + 2*insertSize], readPairs[i][1], reverse=False)
        if minMax3prime != []:

    # Define the anchor point using the minimum-maximum substrings found above

            threePrimeAnchor = np.char.find(refGenome[0][anchorPoint:anchorPoint + 2*insertSize], minMax3prime[0])

    # Look for min-max substrings from the other end (in the opposite direction) upto the anchor point

            minMax5prime = minmaxSubstrings(refGenome[0][anchorPoint + threePrimeAnchor:anchorPoint + threePrimeAnchor + readLength + maxDeletionSize],
                                            readPairs[i][1], reverse=True)
            allMinMaxSubstrings[i].append(1)
            allMinMaxSubstrings[i].append(minMax3prime)
            allMinMaxSubstrings[i].append(minMax5prime)
        else:
            allMinMaxSubstrings[i].append(1)
            allMinMaxSubstrings[i].append([])
            allMinMaxSubstrings[i].append([])

    # If right read is mapped and left is not, look from the left end of right read (in the right-left direction)

    elif np.char.find(refGenome[0], readPairs[i][0]) == -1 and np.char.find(refGenome[0], readPairs[i][1]) > -1:
        anchorPoint = np.char.find(refGenome[0], readPairs[i][1])

        minMax5prime = minmaxSubstrings(refGenome[0][max(0, anchorPoint - 2*insertSize):anchorPoint], readPairs[i][0], reverse=True)
        if minMax5prime != []:

    # Define the anchor point using the minimum-maximum substrings found above

            fivePrimeAnchor = np.char.find(refGenome[0][max(0, anchorPoint - 2*insertSize):anchorPoint][::-1], minMax5prime[0][::-1])

    # Look for min-max substrings from the other end (in the opposite direction) upto the anchor point

            minMax3prime = minmaxSubstrings(refGenome[0][max(0, anchorPoint - fivePrimeAnchor - readLength - maxDeletionSize)
                                                         :anchorPoint - fivePrimeAnchor], readPairs[i][0], reverse=False)
            allMinMaxSubstrings[i].append(0)
            allMinMaxSubstrings[i].append(minMax3prime)
            allMinMaxSubstrings[i].append(minMax5prime)
        else:
            allMinMaxSubstrings[i].append(0)
            allMinMaxSubstrings[i].append([])
            allMinMaxSubstrings[i].append([])

    # If both reads are mapped or none is, there is nothing to be done (except waiting for Godot)

    elif np.char.find(refGenome[0], readPairs[i][1]) > -1 and np.char.find(refGenome[0], readPairs[i][0]) > -1:
        allMinMaxSubstrings[i] == None
    elif np.char.find(refGenome[0], readPairs[i][1]) == -1 and np.char.find(refGenome[0], readPairs[i][0]) == -1:
        allMinMaxSubstrings[i] == None

# Storing the tuples which re-create the original unmapped read

breakPointsDeletions = []
for i in range(totalReadPairs):
    if allMinMaxSubstrings[i] != []:
        if allMinMaxSubstrings[i][1] != [] and allMinMaxSubstrings[i][2] != []:
            if len(allMinMaxSubstrings[i][1][1]) + len(allMinMaxSubstrings[i][2][1]) >= readLength:
                pairIndex = allMinMaxSubstrings[i][0]
                leftBreakPoint = np.char.find(refGenome[0], allMinMaxSubstrings[i][1][1]) + len(allMinMaxSubstrings[i][1][1])
                rightBreakPoint = np.char.find(refGenome[0], allMinMaxSubstrings[i][2][1]) + 0
                breakPointsDeletions.append([leftBreakPoint, rightBreakPoint])

#############################################################
###### Looking for insertions in the reference genome #######
#############################################################

allMinMaxSubstrings = []
for i in range(totalReadPairs):
    allMinMaxSubstrings.append([])

    # If left read is mapped and right is not, start looking from the right end of left read (in the left-right direction)

    if np.char.find(refGenome[0], readPairs[i][0]) > -1 and np.char.find(refGenome[0], readPairs[i][1]) == -1:
        anchorPoint = np.char.find(refGenome[0], readPairs[i][0]) + readLength
        minMax3prime = minmaxSubstrings(refGenome[0][anchorPoint:anchorPoint + 2*insertSize], readPairs[i][1], reverse=False)
        if minMax3prime != []:

    # Define the anchor point using the minimum-maximum substrings found above

            threePrimeAnchor = np.char.find(refGenome[0][anchorPoint:anchorPoint + 2*insertSize], minMax3prime[0])

    # Look for min-max substrings from the other end (in the opposite direction) upto the anchor point

            minMax5prime = minmaxSubstrings(refGenome[0][anchorPoint + threePrimeAnchor:anchorPoint + threePrimeAnchor + readLength - 1],
                                            readPairs[i][1], reverse=True)
            allMinMaxSubstrings[i].append(1)
            allMinMaxSubstrings[i].append(minMax3prime)
            allMinMaxSubstrings[i].append(minMax5prime)
        else:
            allMinMaxSubstrings[i].append(1)
            allMinMaxSubstrings[i].append([])
            allMinMaxSubstrings[i].append([])

    # If right read is mapped and left is not, look from the left end of right read (in the right-left direction)

    elif np.char.find(refGenome[0], readPairs[i][0]) == -1 and np.char.find(refGenome[0], readPairs[i][1]) > -1:
        anchorPoint = np.char.find(refGenome[0], readPairs[i][1])
        minMax5prime = minmaxSubstrings(refGenome[0][max(0, anchorPoint - 2*insertSize):anchorPoint], readPairs[i][0], reverse=True)
        if minMax5prime != []:

    # Define the anchor point using the minimum-maximum substrings found above

            fivePrimeAnchor = np.char.find(refGenome[0][max(0, anchorPoint - 2*insertSize):anchorPoint][::-1], minMax5prime[0][::-1])

    # Look for min-max substrings from the other end (in the opposite direction) upto the anchor point

            minMax3prime = minmaxSubstrings(refGenome[0][max(0, anchorPoint - fivePrimeAnchor - readLength + 1)
                                                         :anchorPoint - fivePrimeAnchor], readPairs[i][0], reverse=False)
            allMinMaxSubstrings[i].append(0)
            allMinMaxSubstrings[i].append(minMax3prime)
            allMinMaxSubstrings[i].append(minMax5prime)
        else:
            allMinMaxSubstrings[i].append(0)
            allMinMaxSubstrings[i].append([])
            allMinMaxSubstrings[i].append([])

    # If both reads are mapped or none is, there is nothing to be done (except waiting for Godot)

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

