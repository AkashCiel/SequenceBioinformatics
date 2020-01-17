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
from NJFunctions import *

"""
This section contains two functions: 
    1. nearestNeighbourExchange 
        As suggested by the name, it performs nearest neighbour interchange for log(number of sequences) times. 
        It re-computes the profiles for the hypothetical joins it considers and uses these profiles to calculate the 
        log-corrected profile distances for each hypothetical join. The join with the minimum log-corrected profile 
        distance is implemented in the tree.
    2. computeBranchLengths
        After nearest neighbour exchange is completed, this function computes the branch length of every branch in the 
        tree using log corrected profile distances.  
"""

def nearestNeighbourExchange(aTree, initialSequences):
    totalRounds = math.log2(initialSequences) + 1

    # Looping until we exhaust total number of NNI rounds
    while totalRounds > 0:
        for node in aTree:

        # Look for an internal node to perform NNI, locate all nodes needed to perform NNI
            if aTree[node][0].getNodeLevel() == "Internal Node":
                fixedNodes = []
                floatingNodes = []
                nodeProfiles = []
                leftNodeChildren = aTree[node][0].getChildren()
                leftNodeParent = aTree[node][0].getParent()
                if leftNodeParent == None:
                    totalRounds = 0
                    break
                rightNodeChildren = aTree[leftNodeParent][0].getChildren()
                rightNodeParent = aTree[leftNodeParent][0].getParent()
                if rightNodeParent == None:
                    totalRounds = 0
                    break
                fixedNodes.append(node)
                fixedNodes.append(leftNodeParent)
                for node in leftNodeChildren:
                    floatingNodes.append(node)
                    nodeProfiles.append(aTree[node][0].getProfile())
                for node in rightNodeChildren:
                    if node not in fixedNodes:
                        floatingNodes.append(node)
                        nodeProfiles.append(aTree[node][0].getProfile())
                floatingNodes.append(rightNodeParent)
                nodeProfiles.append(aTree[rightNodeParent][0].getProfile())

                joinDistances = []

        # Compute all alternative join distances to select the best join from
                joinDistances.append(((-3/4)*(math.log(abs(1 - 4*(profileDelta(nodeProfiles[0], nodeProfiles[1]))/3)))) +
                                     ((-3/4)*(math.log(abs(1 - 4*(profileDelta(nodeProfiles[2], nodeProfiles[3]))/3)))))
                joinDistances.append(((-3/4)*(math.log(abs(1 - 4*(profileDelta(nodeProfiles[0], nodeProfiles[2]))/3)))) +
                                     ((-3/4)*(math.log(abs(1 - 4*(profileDelta(nodeProfiles[1], nodeProfiles[3]))/3)))))
                joinDistances.append(((-3/4)*(math.log(abs(1 - 4*(profileDelta(nodeProfiles[0], nodeProfiles[3]))/3)))) +
                                     ((-3/4)*(math.log(abs(1 - 4*(profileDelta(nodeProfiles[1], nodeProfiles[2]))/3)))))

        # Performing (or not) the neighbour exchange, based on the values of join distances
                if joinDistances.index(min(joinDistances)) == 1:
                    aTree[fixedNodes[0]][0].setChildren([floatingNodes[0], floatingNodes[2]])
                    aTree[fixedNodes[0]][0].setProfile(getAvgProfile(nodeProfiles[0], nodeProfiles[2]))
                    aTree[fixedNodes[1]][0].setChildren([fixedNodes[0], floatingNodes[1]])
                    aTree[fixedNodes[1]][0].setProfile(getAvgProfile(aTree[fixedNodes[0]][0].getProfile(), nodeProfiles[1]))

                    aTree[floatingNodes[1]][0].setParent(fixedNodes[1])
                    aTree[floatingNodes[2]][0].setParent(fixedNodes[0])

                    aTree[fixedNodes[0]][0].updateNNIStatus()
                    aTree[fixedNodes[1]][0].updateNNIStatus()

                elif joinDistances.index(min(joinDistances)) == 2:
                    aTree[fixedNodes[0]][0].setChildren([floatingNodes[0], leftNodeParent])
                    aTree[fixedNodes[0]][0].setParent(rightNodeParent)
                    aTree[fixedNodes[0]][0].setProfile(getAvgProfile(nodeProfiles[0], aTree[leftNodeParent][0].getProfile()))

                    aTree[fixedNodes[1]][0].setChildren([floatingNodes[1], floatingNodes[2]])
                    aTree[fixedNodes[1]][0].setParent(fixedNodes[0])
                    aTree[fixedNodes[1]][0].setProfile(getAvgProfile(nodeProfiles[1], nodeProfiles[2]))

                    aTree[floatingNodes[1]][0].setParent(fixedNodes[1])
                    newChildren = aTree[rightNodeParent][0].getChildren()
                    newChildren.remove(fixedNodes[1])
                    newChildren.append(fixedNodes[0])
                    aTree[rightNodeParent][0].setChildren(newChildren)
                    childProfile01 = aTree[newChildren[0]][0].getProfile()
                    childProfile02 = aTree[newChildren[0]][0].getProfile()
                    aTree[rightNodeParent][0].setProfile(getAvgProfile(childProfile01, childProfile02))

                    aTree[fixedNodes[0]][0].updateNNIStatus()
                    aTree[fixedNodes[1]][0].updateNNIStatus()

        # Reduce total available NNI rounds by 1
                totalRounds -= 1
    return aTree

# Compute branch lengths after fast tree
def computeBranchLengths(aTree):
    allBranchLengths = []
    for node in aTree:

        # Branch length computation if we encounter a leaf node
        if aTree[node][0].getNodeLevel() == "Leaf Node":
            nodeProfile = aTree[node][0].getProfile()
            sibling = (aTree[aTree[node][0].getParent()][0].getChildren()[0] == node) + (aTree[aTree[node][0].getParent()][0].getChildren()[1] == node)
            siblingProfile = aTree[sibling][0].getProfile()
            grandparent = aTree[aTree[node][0].getParent()][0].getParent()
            grandparentProfile = aTree[grandparent][0].getProfile()
            nodeDelta01 = ((-3/4)*(math.log(abs(1 - 4*(profileDelta(nodeProfile, siblingProfile))/3))))
            nodeDelta02 = ((-3/4)*(math.log(abs(1 - 4*(profileDelta(nodeProfile, grandparentProfile))/3))))
            commonDelta = ((-3/4)*(math.log(abs(1 - 4*(profileDelta(siblingProfile, grandparentProfile))/3))))
            branchLength = abs(nodeDelta01 + nodeDelta02 - commonDelta)/2
            allBranchLengths.append([node, aTree[node][0].getParent(), branchLength])

        # Branch length computation if we encounter an internal node
        else:
            if aTree[node][0].getParent() == None:
                continue
            children = aTree[node][0].getChildren()
            childProfile01 = aTree[children[0]][0].getProfile()
            childProfile02 = aTree[children[1]][0].getProfile()
            sibling = (aTree[aTree[node][0].getParent()][0].getChildren()[0] == node) + (aTree[aTree[node][0].getParent()][0].getChildren()[1] == node)
            siblingProfile = aTree[sibling][0].getProfile()

            childSiblingDelta01 = ((-3 / 4) * (math.log(abs(1 - 4 * (profileDelta(childProfile01, siblingProfile)) / 3))))
            childSiblingDelta02 = ((-3 / 4) * (math.log(abs(1 - 4 * (profileDelta(childProfile02, siblingProfile)) / 3))))
            childChildDelta = ((-3 / 4) * (math.log(abs(1 - 4 * (profileDelta(childProfile01, childProfile02)) / 3))))
            grandparent = aTree[aTree[node][0].getParent()][0].getParent()
            if grandparent == None:
                childGrandParentDelta01 = 0.0
                childGrandParentDelta02 = 0.0
                siblingGrandParentDelta = 0.0
            else:
                grandparentProfile = aTree[grandparent][0].getProfile()
                childGrandParentDelta01 = ((-3/4)*(math.log(abs(1 - 4*(profileDelta(childProfile01, grandparentProfile))/3))))
                childGrandParentDelta02 = ((-3/4)*(math.log(abs(1 - 4*(profileDelta(childProfile02, grandparentProfile))/3))))
                siblingGrandParentDelta = ((-3/4)*(math.log(abs(1 - 4*(profileDelta(siblingProfile, grandparentProfile))/3))))

            branchLength = abs((childSiblingDelta01 + childSiblingDelta02 + childGrandParentDelta01 + childGrandParentDelta02)/4 - \
                           (childChildDelta + siblingGrandParentDelta))/2
            allBranchLengths.append([node, aTree[node][0].getParent(), branchLength])
    return allBranchLengths



