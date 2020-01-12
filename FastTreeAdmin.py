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

"""

fastTreeAdmin: FAST TREE ADMINISTRATOR
This function administers the entire fast tree process. There are various functions and classes defined in this project,
each of which perform a smaller operation (for example selecting the best join or re-computing the top hits of a node). 
The fastTreeAdmin function described here orchestrates the entire process of fast tree by calling these functions and 
passing the required parameters as and when needed. 
For the purpose of simplicity, the procedure performed by this function can be seen as comprising of three stages:
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
        b.  Select the best join candidates for this iteration.
        c.  Initialise a parent node for the best join candidates, set them as children of the new node and add this sub-tree
            to the original tree; set age of new node as {age of its older child + 1}
        d.  Inactivate the two children nodes
        e.  Re-compute sum of all up distances
        f.  Re-configure top hits of every node which had one/both of the nodes joined in step b-c
        g.  Compute top hits of the newly added parent node in step c

3. NNI (Nearest Neighbour Interchange):
    During this stage, fastTreeAdmin performs nearest neighbour interchange on the tree constructed in the ITERATIONS stage

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

    theTree = initTree(Sequences)   # Generate ground level nodes
    activeNodes = len(theTree) # Initialise total active nodes count
    totalProfile = computeProfile(Sequences)    # Generate total profile
    iterationNumber = 1 # Keep track of iteration number
    totalUpDistance = allUpDistances(theTree) # Store sum of all up-distances
    for nodeIndex in theTree:    # Generate top hits for each node
        if theTree[nodeIndex][0].getTopHits() == []:
            initialiseTopHits(theTree, theTree[nodeIndex][0], activeNodes, totalProfile, totalUpDistance)

    ##################
    # ITERATIONS STAGE
    ##################

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
        # Set age of new node
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

    ###########
    # NNI STAGE
    ###########

    NNITree = nearestNeighbourExchange(theTree, len(Sequences))
    return NNITree

# Read DNA sequences from text
# Import text file in Python
mainFile: TextIO = open("TestSequences.txt")

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

theTree = fastTreeAdmin(Sequences)
for node in theTree:
    print(theTree[node][0].getIndex(), theTree[node][0].getChildren(), theTree[node][0].getParent())
    """
    if theTree[node][0].getIndex() >= len(Sequences):
        print(theTree[node][0].getIndex(), sequenceID[theTree[node][0].getChildren()[0]],
              sequenceID[theTree[node][0].getChildren()[1]], theTree[node][0].getParent())
    """

"""
Fast tree output:

((ERR2679278:0.00646,(SRR7069540:0.01168,SRR7069811:0.01679)0.950:0.00747)0.972:0.00781,
    ((ERR502513:0.00742,(ERR502505:0.00417,ERR751299:0.00055)0.918:0.00410)0.979:0.01048,
        ((ERR751754:0.00953,(SRR2100230:0.00531,ERR181778:0.00306)0.972:0.00639)0.880:0.00550,
            ((ERR126620:0.00055,
                (ERR403342:0.00916,
                    (SRR6046603:0.00721,
                        ((ERR2307717:0.00215,
                            (SRR1723424:0.03109,ERR552103:0.00286)0.536:0.00129)0.731:0.00325,
                                (SRR6896212:0.00624,(ERR553176:0.00690,ERR369756:0.03306)
                                    0.775:0.00180)
                                        1.000:0.30820)
                                            0.952:0.01401)
                                                0.939:0.00636)
                                                    0.776:0.00119)1.000:0.34104,
            ((ERR181435:0.00055,ERR400415:0.00627)0.826:0.00209,
                ERR234681:0.00312)0.697:0.00053)1.000:0.11717)0.656:0.00716)0.991:0.01189,
    ((ERR2513433:0.01121,SRR5657468:0.00144)1.000:0.03326,ERR017782:0.04824)0.340:0.00359);
    
"""
