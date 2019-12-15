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

# Read DNA sequences from text
# Import text files in Python
reads01: TextIO = open("PindelReads01.txt")
reads02: TextIO = open("PindelReads02.txt")
referenceGenome: TextIO = open("PindelReferenceGenome.txt")

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

for i in range(50):
    print(i, np.char.find(refGenome[0], leftEndReads[i]), np.char.find(refGenome[0], rightEndReads[i]))
    if np.char.find(refGenome[0], leftEndReads[i]) == -1:
        print(leftEndReads[i])
        print(minmaxSubstrings(refGenome[0], leftEndReads[i], reverse=False))
        print(minmaxSubstrings(refGenome[0], leftEndReads[i], reverse=True))
        print("Reversed")
        print(minmaxSubstrings(refGenome[0], leftEndReads[i][::-1], reverse=False))
        print(minmaxSubstrings(refGenome[0], leftEndReads[i][::-1], reverse=True))
    if np.char.find(refGenome[0], rightEndReads[i]) == -1:
        print(rightEndReads[i])
        print(minmaxSubstrings(refGenome[0], rightEndReads[i], reverse=False))
        print(minmaxSubstrings(refGenome[0], rightEndReads[i], reverse=True))
        print("Reversed")
        print(minmaxSubstrings(refGenome[0], rightEndReads[i][::-1], reverse=False))
        print(minmaxSubstrings(refGenome[0], rightEndReads[i][::-1], reverse=True))