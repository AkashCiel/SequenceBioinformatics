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

# Read DNA sequences from text
# Import text file in Python

##################################################################
################# Type name of input file below ##################
##################################################################
mainFile: TextIO = open("4Sequences.txt")

"""

fastTreeAdmin: FAST TREE ADMINISTRATOR
This function administers the entire fast tree process. There are various functions and classes defined in this project,
each of which perform a smaller operation (for example selecting the best join or re-computing the top hits of a node). 
The fastTreeAdmin function described here orchestrates the entire process of fast tree by calling these functions and 
passing the required parameters as and when needed. 
For the purpose of simplicity, the procedure performed by this function can be seen as comprising of four stages:
1. INITIALISATION: 
    Operations performed during this stage are:
        a.  Initialise all the nodes of the tree and assign a profile from the input sequences to each (more details on 
            'node' class in FastTreeClasses documentation
        b.  Generate the total profile of input sequences, used during the computation of neighbour-joining criterion
        c.  Calculate the initial total up distance of the tree, to be used for the first join
        d.  Generate top hits for each of the initial nodes       
        
2. ITERATIONS:
    During the iterations stage, fastTreeAdmin begins joining nodes and building the tree in an iterative fashion. 
    During this stage, the following operations are performed iteratively until we are left with only 2 active nodes (a node
    which has been joined with another node is inactivated and is not considered for future joins). Operations 
    performed during this stage are:
        a.  Increment the age of all active nodes by 1. This ensures that once a node is too old, its top hits are 
            re-computed. 
        b.  Select the best join candidates from the list of top-hits for each node for this iteration (more details about
            finding the best join in NJFunctions documentation).
        c.  Initialise a parent node for the best join candidates, set them as children of the new node and add this sub-tree
            to the original tree; set age of new node as {age of its older child + 1}
        d.  Inactivate the two children nodes
        e.  Re-compute sum of all up distances
        f.  Re-configure top hits of every node which had one/both of the nodes joined in step b-c
        g.  Compute top hits of the newly added parent node in step c

3. NNI (Nearest Neighbour Interchange):
    During this stage, fastTreeAdmin performs nearest neighbour interchange on the tree constructed in the ITERATIONS stage
    
4. BRANCH LENGTH Computation:
    During this stage, fastTreeAdmin uses the computeBranchLength function to calculate the log-corrected branch lengths
    for the entire constructed tree

PLEASE NOTE: 
    All the computations mentioned above are performed by some function defined specifically for that task. These 
    functions can be found in the files NJFunctions.py and NNIFunctions.py, along with function specific documentation. 
    The code in this section, as mentioned above, simply orchestrates the entire process of fast tree

"""
# FastTree Administrator
def fastTreeAdmin(Sequences):

    ######################
    # INITIALISATION STAGE
    ######################

    # Generate ground level nodes
    theTree = initTree(Sequences)

    # Initialise total active nodes count
    activeNodes = len(theTree)

    # Generate total profile
    totalProfile = computeProfile(Sequences)

    # Keep track of iteration number
    iterationNumber = 1

    # Store sum of all up-distances
    totalUpDistance = allUpDistances(theTree)

    # Generate top hits for each node which does not yet has a top-hits list
    for nodeIndex in theTree:
        if theTree[nodeIndex][0].getTopHits() == []:
            initialiseTopHits(theTree, theTree[nodeIndex][0], activeNodes, totalProfile, totalUpDistance)

    ##################
    # ITERATIONS STAGE
    ##################

    while activeNodes > 2:

        # Increment all active node ages by 1
        for node in theTree:
            if theTree[node][0].getNodeStatus() == "Active":
                nodeAge = theTree[node][0].getAge()
                theTree[node][0].setAge(nodeAge + 1)

        # Find the best join for this iteration
        bestJoinTuple = getBestJoin(theTree, activeNodes, totalProfile, totalUpDistance)[:2]

        # Initialise a new node with the appropriate profile and node index
        newNodeIndex = len(theTree)
        newNode = fastTreeNode(newNodeIndex, getAvgProfile(theTree[bestJoinTuple[0]][0].getProfile(), theTree[bestJoinTuple[1]][0].getProfile()))

        # Add new node in the tree
        theTree[newNodeIndex].append(newNode)
        theTree[newNodeIndex][0].setChildren(bestJoinTuple)

        # Set age of new node
        theTree[newNodeIndex][0].setAge(max(theTree[bestJoinTuple[0]][0].getAge(), theTree[bestJoinTuple[1]][0].getAge()) + 1)

        # Establish child-parent connections among the proper nodes and deactivate the children nodes
        theTree[bestJoinTuple[0]][0].setParent(newNodeIndex)
        theTree[bestJoinTuple[0]][0].inActivate()
        theTree[bestJoinTuple[1]][0].setParent(newNodeIndex)
        theTree[bestJoinTuple[1]][0].inActivate()

        # Reduce active node count by 1
        activeNodes -= 1

        # Re-compute sum of all up-distances
        if activeNodes > 2:
            totalUpDistance = allUpDistances(theTree)
            reconfigureTopHits(theTree, bestJoinTuple, newNodeIndex, activeNodes, totalProfile, totalUpDistance)

    # Look for the last join if we are left with only 2 active nodes
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

    ###########
    # NNI STAGE
    ###########

    NNITree = nearestNeighbourExchange(theTree, len(Sequences))
    return NNITree


# Read file into a text string
InputText = mainFile.readlines()
Sequences = []
sequenceID = []
for i in range(0, len(InputText)):
    InputText[i] = InputText[i].rstrip('\n')
    if (i % 2) == 0:
        sequenceID.append(InputText[i])
    else:
        Sequences.append(InputText[i])

# Generate a phylogenetic tree post nearest-neighbour interchange
theTree = fastTreeAdmin(Sequences)

#########################################################
############# BRANCH LENGTH COMPUTATION #################
#########################################################

branchLengths = computeBranchLengths(theTree)

###############################################################
######### PRINTING THE TREE NODES AND BRANCH LENGTHS ##########
###############################################################

for node in theTree:
    if theTree[node][0].getNodeLevel() == "Leaf Node":
        print("Index: ", theTree[node][0].getIndex(), ", Input ID: ", sequenceID[theTree[node][0].getIndex()], ", No children nodes")
        print("Parent Index: ", theTree[node][0].getParent(), ", Branch Length = ", branchLengths[theTree[node][0].getIndex()][2])
    else:
        print("Index: ", theTree[node][0].getIndex())
        print(" Child 1 Index: ", theTree[node][0].getChildren()[0], ", Branch Length = ", branchLengths[theTree[node][0].getChildren()[0]][2])
        print(" Child 2 Index: ", theTree[node][0].getChildren()[1], ", Branch Length = ", branchLengths[theTree[node][0].getChildren()[1]][2])
        if theTree[node][0].getParent() != None:
            print("Parent Index: ", theTree[node][0].getParent(), ", Branch Length = ", branchLengths[theTree[node][0].getIndex()][2])
        else:
            print("No parent")
