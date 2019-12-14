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
        self.parent = None
        self.leftChild = None
        self.rightChild = None
        self.age = 0
    def inActivate(self):
        self.status = "Inactive"
    def getAge(self):
        return self.age
    def setAge(self, newAge):
        self.age = newAge
    def getNodeLevel(self):
        if self.leftChild == None and self.rightChild == None: return "Leaf Node"
        else: return "Non-leaf Node"
    def getProfile(self):
        return self.profile
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
    def setLeftChild(self, index):
        self.leftChild = index
    def setRightChild(self, index):
        self.rightChild = index
    def setParent(self, index):
        self.parent = index
    def getRelatives(self):
        return [self.parent, self.leftChild, self.rightChild]
