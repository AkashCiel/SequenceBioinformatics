from typing import TextIO
from operator import itemgetter
from random import seed
from random import randint
from collections import defaultdict
from copy import deepcopy
from random import randint
from decimal import *
import math

from FastTreeClasses import *
"""
This section contains all the functions needed for the INITIALISATION and ITERATIONS stages of fast tree. Description is as follows:

    1. patternToNumber
        Converts a pattern into its lexicographical order number
    2. computeProfile
        Computes the profile of a given set of DNA sequences and returns a multi-dimensional array containing these profiles
    3. profileDelta
        Calculates the sum of differences in probabilities at all positions of two given profiles. (There was some ambiguity
        about whether the sum of differences should be divided by the length of sequence or not. I tried both and the results
        seem to be 'more correct' if I do not divide by the length of sequence. Therefore, this implementation does not
        divide the sum of differences by the length of the sequence.
    4. getAvgProfile
        Calculates the average profile based on two input profiles. Used for calculating profile of a parent node given 
        the profiles of its two children
    5. initTree
        Initialises a tree given a set of DNA sequences with each sequence represented by an instance of the class fastTreeNode
        (more details about this class in FastTreeClasses documentation)
    6. getUpDistance
        Calculates the up-distance of a given node in a given tree
    7. allUpDistances
        Calculates the sum of up-distances of all nodes in a given tree
    8. neighbourJoinScore
        Calculates the neigbour joining score of two given nodes in a given tree. Uses the function getUpDistance to calculate
        individual up-distances
    9. initialiseTopHits
        Initialises the top hits for a given node. Additionally, it then intialises the top-hits for the top 'm' hits of the 
        original node where 'm' is equal to square root of total sequences in the tree
    10. refreshTopHits
        This function is called if a node has too few top-hits (because many nodes from its top-hits have been joined and
        hence, inactivated) or if a node is too 'old' (determined by a node's incrementing age as explained in FastTreeAdmin
        documentation). In either case, refreshTopHits 'replenishes' the list of top-hits for a given node.
    11. getBestJoin
        This function is called in each iteration of fast tree to select a best join pair. It first selects top m best-joins
        from the top-hits list of each active node. Then, it updates these top m best-joins by recomputing the neighbour 
        joining score. It then picks the best join (A, B) from this updated list and then performs a local hill climbing search 
        by considering the join of node A with all top-hits of node B and vice versa. The best join obtained through this 
        elaborate process is returned to fastTreeAdmin to be joined in a given iteration.
    12. getBestJoinV02
        This function is called when there are two active nodes left. It simply returns those two active nodes to be joined.
    13. reconfigureTopHits
        This function is called every time a join is performed. It replaces the two recently joined nodes from the top-hits
        list of every node with the new parent node and generates a top-hits list of the parent node using the top-hits
        lists of its two (now inactive) children. Finally, it checks if the parent node has too few top-hits or is too old.
        In either case, it calls refreshTopHits to 'replenish' the list of top-hits of the parent node.  
"""
# Patter to number
def patternToNumber(pattern):
    k = len(pattern)
    if k == 1:
        if pattern == 'A': return 0
        elif pattern == 'C': return 1
        elif pattern == 'G': return 2
        elif pattern == 'T': return 3
    else:
        return 4*patternToNumber(pattern[:k-1]) + patternToNumber(pattern[k-1])

# Compute profile
def computeProfile(Sequences):
    TotalProfile = []
    totalSequences = len(Sequences)

    #Initialise Total Profile matrix
    for i in range(len(Sequences[0])):
        TotalProfile.append([0,0,0,0])

    #Fill in the profile matrix
    for aSequence in Sequences:
        for i in range(len(aSequence)):
            if aSequence[i] != '-':
                TotalProfile[i][patternToNumber(aSequence[i])] += (1/totalSequences)
    return TotalProfile

# Compute profile delta
def profileDelta(profile01, profile02):
    totalDelta = 0
    totalNucleotides = len(profile01)
    for i in range(totalNucleotides):
        for j in range(4):
            totalDelta += abs(profile01[i][j] - profile02[i][j])
    return totalDelta

# Compute Average Profile
def getAvgProfile(profile01, profile02):
    avgProfile = []
    for vectorIndex in range(len(profile01)):
        avgProfile.append([0,0,0,0])
        for i in range(4):
            avgProfile[vectorIndex][i] = (profile01[vectorIndex][i] + profile02[vectorIndex][i])/2
    return avgProfile

# Initialise tree
def initTree(Sequences):
    nodeNumber = 0
    theTree = defaultdict(list)
    for aSequence in Sequences:
        aNode = fastTreeNode(nodeNumber, computeProfile([aSequence]))
        theTree[nodeNumber].append(aNode)
        nodeNumber += 1
    return theTree

# Get up distance
def getUpDistance(aTree, aNodeIndex):
    if aTree[aNodeIndex][0].getNodeLevel() == "Leaf Node":
        return 0
    else:
        leftChild = aTree[aNodeIndex][0].getChildren()[0]
        rightChild = aTree[aNodeIndex][0].getChildren()[1]
        return ((profileDelta(aTree[leftChild][0].getProfile(), aTree[rightChild][0].getProfile())) / 2)

# Store all up-distances
def allUpDistances(aTree):
    sumUpDistances = 0
    for node in aTree:
        if aTree[node][0].getNodeStatus() == "Active":
            sumUpDistances += getUpDistance(aTree, node)
    return sumUpDistances

# Compute Neighbour Join Score
def neighbourJoinScore(aTree, indexNode01, indexNode02, totalActiveNodes, totalProfile, totalUpDist):
    profileNode01 = aTree[indexNode01][0].getProfile()
    profileNode02 = aTree[indexNode02][0].getProfile()

    # Compute up-distances and out-distances for both nodes
    updist01 = getUpDistance(aTree, indexNode01)
    updist02 = getUpDistance(aTree, indexNode02)
    outdist01 = ((totalActiveNodes*profileDelta(profileNode01, totalProfile)) - (2*updist01) - (totalActiveNodes-2)*updist01 - totalUpDist)/(totalActiveNodes - 2)
    outdist02 = ((totalActiveNodes * profileDelta(profileNode02, totalProfile)) - (2*updist02) - (totalActiveNodes-2)*updist02 - totalUpDist) / (totalActiveNodes - 2)

    # Compute the join score based on up-distances, out-distances, and profile difference
    joinScore = profileDelta(profileNode01, profileNode02) - updist01 - updist02 - outdist01 - outdist02
    return joinScore

# Initialise top hits for a node and its neighbours
def initialiseTopHits(aTree, aNode, totalActiveNodes, totalProfile, totalUpDist):
    totalSequences = len(aTree)
    totalHits = 2*int(math.sqrt(totalSequences))
    topHits = []

    # Compute top 2m hits for the given node
    while len(topHits) <= (totalHits):
        for i in range(totalSequences):
            if i != aNode.getIndex():
                topHits.append([i, int(neighbourJoinScore(aTree, aNode.getIndex(), i, totalActiveNodes, totalProfile, totalUpDist))])
    topHits.sort(key=itemgetter(1))
    topHits = topHits[:totalHits]
    aNode.addTopHits(topHits, totalHits)

    # Add top hits for top m hits of the given node
    # Consider top hits of the original node first
    for i in range(int(totalHits/2)):
        topHits = []
        selfHit = [aNode.getIndex(), aNode.getTopHits()[i][1]]
        topHits.append(selfHit)
        neighbourIndex = aNode.getTopHits()[i][0]
        for aHit in aNode.getTopHits():
            if aHit[0] != neighbourIndex:
                topHits.append([aHit[0], int(neighbourJoinScore(aTree, neighbourIndex, aHit[0], totalActiveNodes, totalProfile, totalUpDist))])
        topHits.sort(key=itemgetter(1))
        topHits = topHits[:totalHits]
        aTree[neighbourIndex][0].addTopHits(topHits, totalHits)

# Re-plenish Top Hits for a node with too few top hits or for a node which is too old
def refreshTopHits(aTree, aNode, totalActiveNodes, totalProfile, totalUpDist):
    totalSequences = len(aTree)
    totalHits = 2*int(math.sqrt(totalSequences))
    topHits = []
    for i in range(totalSequences):
        if i != aTree[aNode][0].getIndex() and aTree[i][0].getNodeStatus() == "Active":
            topHits.append([i, int(neighbourJoinScore(aTree, aTree[aNode][0].getIndex(), i, totalActiveNodes, totalProfile, totalUpDist))])
    topHits.sort(key=itemgetter(1))
    topHits = topHits[:totalHits]
    aTree[aNode][0].addTopHits(topHits, totalHits)
    aTree[aNode][0].setAge(0)
    newNodeTopHits = aTree[aNode][0].getTopHits()

    #Compare all top hit nodes against each other
    allCombinationHits = defaultdict(list)
    for i in range(len(newNodeTopHits)):
        for j in range(len(newNodeTopHits)):
            node1 = newNodeTopHits[i][0]
            node2 = newNodeTopHits[j][0]
            if i != j and any(node2 in sublist for sublist in allCombinationHits[node1]) == False:
                NJscore = int(neighbourJoinScore(aTree, node1, node2, totalActiveNodes, totalProfile, totalUpDist))
                allCombinationHits[node1].append([node2, NJscore])
                allCombinationHits[node2].append([node1, NJscore])

    #Add top hits to each node
    for node in allCombinationHits:
        aTree[node][0].addTopHits(allCombinationHits[node], totalHits)
        aTree[node][0].setAge(0)

# Search for best join in a tree
def getBestJoin(aTree, totalActiveNodes, totalProfile, totalUpDist):
    totalJoins = int(math.sqrt(totalActiveNodes))
    mBestJoins = []

    # Find m best joins
    for node in aTree:
        if aTree[node][0].getNodeStatus() == "Active":
            if [aTree[node][0].getTopHits()[0][0], node, aTree[node][0].getTopHits()[0][1]] not in mBestJoins:
                mBestJoins.append([node, aTree[node][0].getTopHits()[0][0], aTree[node][0].getTopHits()[0][1]])
    mBestJoins.sort(key=itemgetter(2))
    mBestJoins = mBestJoins[:totalJoins]

    # Update neighbour join score of m best joins
    for i in range(len(mBestJoins)):
        mBestJoins[i][2] = neighbourJoinScore(aTree, mBestJoins[i][0], mBestJoins[i][1], totalActiveNodes, totalProfile, totalUpDist)

    # Sort by updated neighbour score
    mBestJoins.sort(key=itemgetter(2))
    mBestJoins = mBestJoins[0]

    # Perform Local hill climbing search for an even better best join
    bestJoinTuple = [mBestJoins]
    for i in range(2):
        for node in aTree[mBestJoins[i]][0].getTopHits():
            neighbourJoinDistance = neighbourJoinScore(aTree, mBestJoins[i], node[0], totalActiveNodes, totalProfile, totalUpDist)
            bestJoinTuple.append([mBestJoins[i], node[0], neighbourJoinDistance])
    bestJoinTuple.sort(key=itemgetter(2))

    # Return the best join tuple
    return bestJoinTuple[0]

# Search for last 2 active nodes
def getBestJoinV02(aTree):
    bestJoinTuple = []
    for node in aTree:
        if aTree[node][0].getNodeStatus() == "Active":
            bestJoinTuple.append(aTree[node][0].getIndex())
    return bestJoinTuple

# Reconfigure top hits after a join
def reconfigureTopHits(theTree, bestJoinTuple, newNodeIndex, totalActiveNodes, totalProfile, totalUpDist):
    totalHits = 2 * int(math.sqrt(totalActiveNodes))

    # Replace the now inactive children by their ancestor node
    for node in theTree:

        # Remove the recently joined nodes from top-hits list of nodes, wherever they appear
        if theTree[node][0].getIndex() not in bestJoinTuple:
            nodeTopHits = theTree[node][0].getTopHits()
            newTopHits = []
            for sublist in nodeTopHits:
                if sublist[0] not in bestJoinTuple:
                    newTopHits.append(sublist)
                else:
                    newTopHits.append([newNodeIndex, neighbourJoinScore(theTree, theTree[node][0].getIndex(), newNodeIndex, totalActiveNodes, totalProfile, totalUpDist)])
            theTree[node][0].clearTopHits()
            theTree[node][0].addTopHits(newTopHits, totalHits)

        # Reset age of the updated node to 0
            theTree[node][0].setAge(0)

    # Generate top-hits for the parent using top-hits of its children
    parentTopHits = []
    leftChildHits = theTree[bestJoinTuple[0]][0].getTopHits()
    for sublist in leftChildHits:
        joinScore = neighbourJoinScore(theTree, newNodeIndex, sublist[0], totalActiveNodes, totalProfile, totalUpDist)
        parentTopHits.append([sublist[0], joinScore])
    rightChildHits = theTree[bestJoinTuple[1]][0].getTopHits()
    for sublist in rightChildHits:
        joinScore = neighbourJoinScore(theTree, newNodeIndex, sublist[0], totalActiveNodes, totalProfile, totalUpDist)
        parentTopHits.append([sublist[0], joinScore])

    correctedHits = []
    for sublist in parentTopHits:
        if theTree[sublist[0]][0].getNodeStatus() == "Active": correctedHits.append(sublist)
    theTree[newNodeIndex][0].addTopHits(correctedHits, totalHits)

    # Refresh top hits of parent node if it has too few hits or is too old
    if len(theTree[newNodeIndex][0].getTopHits()) < 0.4*totalHits or theTree[newNodeIndex][0].getAge() > (math.log2(totalHits)):
        refreshTopHits(theTree, newNodeIndex, totalActiveNodes, totalProfile, totalUpDist)
