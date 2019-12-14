from typing import TextIO
from operator import itemgetter
from random import seed
from random import randint
from collections import defaultdict
from copy import deepcopy
from random import randint
from decimal import *
import math

from NJFunctions import *
from FastTreeClasses import *
from NNIFunctions import *

# FastTree Administrator
def fastTreeAdmin(Sequences):
    theTree = initTree(Sequences)   # Generate ground level nodes
    activeNodes = len(theTree) # Initialise total active nodes count
    totalProfile = computeProfile(Sequences)    # Generate total profile
    iterationNumber = 1 # Keep track of iteration number
    totalUpDistance = allUpDistances(theTree) # Store sum of all up-distances
    for nodeIndex in theTree:    # Generate top hits for each node
        if theTree[nodeIndex][0].getTopHits() == []:
            initialiseTopHits(theTree, theTree[nodeIndex][0], activeNodes, totalProfile, totalUpDistance)
    while activeNodes > 2:
        for node in theTree: # Increment all active node ages by 1
            if theTree[node][0].getNodeStatus() == "Active":
                nodeAge = theTree[node][0].getAge()
                theTree[node][0].setAge(nodeAge + 1)
        bestJoinTuple = getBestJoin(theTree, activeNodes, totalProfile, totalUpDistance)[:2] # Find the best join
        newNodeIndex = len(theTree)
        newNode = fastTreeNode(newNodeIndex, getAvgProfile(theTree[bestJoinTuple[0]][0].getProfile(), theTree[bestJoinTuple[1]][0].getProfile()))
        theTree[newNodeIndex].append(newNode) # Add new node in the tree
        theTree[newNodeIndex][0].setChildren(bestJoinTuple)
        theTree[newNodeIndex][0].setAge(max(theTree[bestJoinTuple[0]][0].getAge(), theTree[bestJoinTuple[1]][0].getAge()) + 1)
        theTree[bestJoinTuple[0]][0].setParent(newNodeIndex)
        theTree[bestJoinTuple[0]][0].inActivate()
        theTree[bestJoinTuple[1]][0].setParent(newNodeIndex)
        theTree[bestJoinTuple[1]][0].inActivate()
        activeNodes -= 1
        if activeNodes > 2:
            totalUpDistance = allUpDistances(theTree)  # Re-compute sum of all up-distances
            reconfigureTopHits(theTree, bestJoinTuple, newNodeIndex, activeNodes, totalProfile, totalUpDistance)
    if activeNodes == 2:
        bestJoinTuple = getBestJoinV02(theTree)
        newNodeIndex = len(theTree)
        newNode = fastTreeNode(newNodeIndex, getAvgProfile(theTree[bestJoinTuple[0]][0].getProfile(), theTree[bestJoinTuple[1]][0].getProfile()))
        theTree[newNodeIndex].append(newNode)  # Add new node in the tree
        theTree[newNodeIndex][0].setChildren(bestJoinTuple)
        theTree[bestJoinTuple[0]][0].setParent(newNodeIndex)
        theTree[bestJoinTuple[0]][0].inActivate()
        theTree[bestJoinTuple[1]][0].setParent(newNodeIndex)
        theTree[bestJoinTuple[1]][0].inActivate()
    return theTree

# Read DNA sequences from text
# Import text file in Python
mainFile: TextIO = open("TestSequences.txt")

# Read file into a text string
Sequences = mainFile.readlines()
for i in range(0, len(Sequences)):
    Sequences[i] = Sequences[i].rstrip('\n')

theTree = fastTreeAdmin(Sequences)
for node in theTree:
    print(theTree[node][0].getIndex(), theTree[node][0].getChildren(), theTree[node][0].getParent())

NNITree = nearestNeighbourExchange(theTree, 10)
for node in NNITree:
    print(NNITree[node][0].getIndex(), NNITree[node][0].getChildren(), NNITree[node][0].getParent())
