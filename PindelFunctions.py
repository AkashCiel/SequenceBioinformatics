"""
The Pindel program uses only one function, minmaxSubstrings. As suggested by the name, this function returns a tuple.
This tuple is the minimum and maximum unique substrings of "pattern" as found in the input "genome". The function takes an
additional parameter called "reverse" which takes either "True" or "False" as acceptable values. If the value of "reverse"
is set to "False", the function will consider substrings from the left end of "pattern" and look for minimum and maximum
unique substrings starting from the left end of "genome". If the value of "reverse" is set to be true, the function will
consider substrings from the right end of "pattern" and look for minimum and maximum unique substrings starting from the
right end of "genome". While returning the minimum-maximum tuple though, the substrings should always be read in the
left-to-right direction.
Thus, in case of a deletion event, this program will obtain minimum-maximum tuple of a read for both values of "reverse"
(True and False). Combining these substrings will result in perfect reproduction of the original read.
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

# Find minimum and maximum unique substrings
def minmaxSubstrings(genome, pattern, reverse):

    # Reversing the input strings if "reverse" is set at "True"
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
