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

def nearestNeighbourExchange(aTree, initialSequences):
    totalRounds = math.log2(initialSequences) + 1
    while totalRounds > 0:
        for node in aTree:
            if aTree[node][0].getNodeLevel() == "Internal Node" and aTree[node][0].getNNIStatus() == False:
                fixedNodes = []
                floatingNodes = []
                nodeProfiles = []
                leftNodeChildren = aTree[node][0].getChildren()
                leftNodeParent = aTree[node][0].getParent()
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
                joinDistances.append(((-3/4)*(math.log(abs(1 - 4*(profileDelta(nodeProfiles[0], nodeProfiles[1]))/3)))) +
                                     ((-3/4)*(math.log(abs(1 - 4*(profileDelta(nodeProfiles[2], nodeProfiles[3]))/3)))))
                joinDistances.append(((-3/4)*(math.log(abs(1 - 4*(profileDelta(nodeProfiles[0], nodeProfiles[2]))/3)))) +
                                     ((-3/4)*(math.log(abs(1 - 4*(profileDelta(nodeProfiles[1], nodeProfiles[3]))/3)))))
                joinDistances.append(((-3/4)*(math.log(abs(1 - 4*(profileDelta(nodeProfiles[0], nodeProfiles[3]))/3)))) +
                                     ((-3/4)*(math.log(abs(1 - 4*(profileDelta(nodeProfiles[1], nodeProfiles[2]))/3)))))

                if joinDistances.index(min(joinDistances)) == 1:
                    print(fixedNodes, 1)
                    aTree[fixedNodes[0]][0].setChildren([floatingNodes[0], floatingNodes[2]])
                    aTree[fixedNodes[0]][0].setProfile(getAvgProfile(nodeProfiles[0], nodeProfiles[2]))
                    aTree[fixedNodes[1]][0].setChildren([fixedNodes[0], floatingNodes[1]])
                    aTree[fixedNodes[1]][0].setProfile(getAvgProfile(aTree[fixedNodes[0]][0].getProfile(), nodeProfiles[1]))

                    aTree[floatingNodes[1]][0].setParent(fixedNodes[1])
                    aTree[floatingNodes[2]][0].setParent(fixedNodes[0])

                    aTree[fixedNodes[0]][0].updateNNIStatus()
                    aTree[fixedNodes[1]][0].updateNNIStatus()

                elif joinDistances.index(min(joinDistances)) == 2:
                    print(fixedNodes, 2)
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

                totalRounds -= 1
    return aTree





