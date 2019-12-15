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
def minmaxSubstrings(genome, pattern):
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
    finalSubStrings.append(subStrings[0])
    finalSubStrings.append(subStrings[len(subStrings) - 1])
    return finalSubStrings

