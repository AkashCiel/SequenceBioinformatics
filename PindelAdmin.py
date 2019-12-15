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

refGenome = referenceGenome.readlines()
for i in range(0, len(refGenome)):
    Sequences[i] = Sequences[i].rstrip('\n')

values = np.array([1,2,3,1,2,4,5,6,3,2,1])
searchval = 3
ii = np.where('A' in leftEndReads[0])

print(minmaxSubstrings(leftEndReads[0], 'GAGGAGTCACT'))
#print('ABC' in 'THUINAHABCKMJ')

