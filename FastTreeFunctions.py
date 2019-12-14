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

# Patter to number
def patternToNumber(pattern):
    k = len(pattern)
    if k == 1:
        if pattern == 'A' or pattern == 'N': return 0
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
        leftChild = aTree[aNodeIndex][0].leftChild
        rightChild = aTree[aNodeIndex][0].rightChild
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
    updist01 = getUpDistance(aTree, indexNode01)
    updist02 = getUpDistance(aTree, indexNode02)
    outdist01 = ((totalActiveNodes*profileDelta(profileNode01, totalProfile)) - (totalActiveNodes-2)*updist01 - totalUpDist)/(totalActiveNodes - 2)
    outdist02 = ((totalActiveNodes * profileDelta(profileNode02, totalProfile)) - (totalActiveNodes - 2)*updist02 - totalUpDist) / (totalActiveNodes - 2)
    joinScore = profileDelta(profileNode01, profileNode02) - updist01 - updist02 - outdist01 - outdist02
    return joinScore

# Initialise top hits for a node and its neighbours
def initialiseTopHits(aTree, aNode, totalActiveNodes, totalProfile, totalUpDist):
    totalSequences = len(aTree)
    totalHits = 2*int(math.sqrt(totalSequences))
    topHits = []
    while len(topHits) <= (totalHits):
        for i in range(totalSequences):
            if i != aNode.getIndex():
                topHits.append([i, int(neighbourJoinScore(aTree, aNode.getIndex(), i, totalActiveNodes, totalProfile, totalUpDist))])
    topHits.sort(key=itemgetter(1))
    topHits = topHits[:totalHits]
    aNode.addTopHits(topHits, totalHits)
    # Add top hits for top m hits of aNode
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

# Fix Top Hits for a node with too few top hits
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

# Search for best join in a tree
def getBestJoin(aTree, totalActiveNodes, totalProfile, totalUpDist):
    totalJoins = int(math.sqrt(totalActiveNodes))
    mBestJoins = []
    for node in aTree: # Find m best joins
        if aTree[node][0].getNodeStatus() == "Active":
            if [aTree[node][0].getTopHits()[0][0], node, aTree[node][0].getTopHits()[0][1]] not in mBestJoins:
                mBestJoins.append([node, aTree[node][0].getTopHits()[0][0], aTree[node][0].getTopHits()[0][1]])
    mBestJoins.sort(key=itemgetter(2))
    mBestJoins = mBestJoins[:totalJoins]
    #print(mBestJoins)
    for i in range(len(mBestJoins)): # Update neighbour join score of m best joins
        mBestJoins[i][2] = neighbourJoinScore(aTree, mBestJoins[i][0], mBestJoins[i][1], totalActiveNodes, totalProfile, totalUpDist)
    mBestJoins.sort(key=itemgetter(2)) # Sort by updated neighbour score
    mBestJoins = mBestJoins[0]
    #Local hill climbing search
    bestJoinTuple = [mBestJoins]
    for i in range(2):
        for node in aTree[mBestJoins[i]][0].getTopHits():
            neighbourJoinDistance = neighbourJoinScore(aTree, mBestJoins[i], node[0], totalActiveNodes, totalProfile, totalUpDist)
            bestJoinTuple.append([mBestJoins[i], node[0], neighbourJoinDistance])
    bestJoinTuple.sort(key=itemgetter(2))

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
    for node in theTree: # Replace eliminated children by their ancestor node
        if theTree[node][0].getIndex() not in bestJoinTuple: # Do not update best join tuple, they've been inactivated
            nodeTopHits = theTree[node][0].getTopHits()
            newTopHits = []
            for sublist in nodeTopHits:
                if sublist[0] not in bestJoinTuple:
                    newTopHits.append(sublist)
                else:
                    newTopHits.append([newNodeIndex, neighbourJoinScore(theTree, theTree[node][0].getIndex(), newNodeIndex, totalActiveNodes, totalProfile, totalUpDist)])
            theTree[node][0].clearTopHits()
            theTree[node][0].addTopHits(newTopHits, totalHits)

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
    if len(theTree[newNodeIndex][0].getTopHits()) < totalHits:
        refreshTopHits(theTree, newNodeIndex, totalActiveNodes, totalProfile, totalUpDist)
