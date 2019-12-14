from typing import TextIO
from operator import itemgetter
from random import seed
from random import randint
from collections import defaultdict
from copy import deepcopy
from random import randint
from decimal import *
import math

# Define node class
class fastTreeNode:
    def __init__(self, index, profile):
        self.index = index
        self.profile = profile
        self.status = "Active"
        self.topHits = []
        self.children = []
        self.parent = None
        self.age = 0
        self.nniStatus = False
    def inActivate(self):
        self.status = "Inactive"
    def getAge(self):
        return self.age
    def setAge(self, newAge):
        self.age = newAge
    def getNodeLevel(self):
        if self.children == []: return "Leaf Node"
        else: return "Internal Node"
    def getProfile(self):
        return self.profile
    def setProfile(self, newProfile):
        self.profile = newProfile
    def getNodeStatus(self):
        return self.status
    def getTopHits(self):
        return self.topHits
    def addTopHits(self, newTopHits, totalHits):
        for aHit in newTopHits:
            if any(aHit[0] in sublist for sublist in self.topHits) == False:
                self.topHits.append(aHit)
        self.topHits.sort(key=itemgetter(1))
        self.topHits = self.topHits[:totalHits]
    def clearTopHits(self):
        self.topHits.clear()
    def getIndex(self):
        return self.index
    def setChildren(self, children):
        self.children = children
    def getChildren(self):
        return self.children
    def setParent(self, parent):
        self.parent = parent
    def getParent(self):
        return self.parent
    def getNNIStatus(self):
        return self.nniStatus
    def updateNNIStatus(self):
        self.nniStatus = True
