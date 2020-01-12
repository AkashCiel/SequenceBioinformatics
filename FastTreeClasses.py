from typing import TextIO
from operator import itemgetter
from random import seed
from random import randint
from collections import defaultdict
from copy import deepcopy
from random import randint
from decimal import *
import math

"""
fastTreeNode: a class defined for the purpose of fast tree
    A fast tree node is meant to serve as the data structure storing details of a node in the phylogenic tree constructed
    by fast tree. Variables of this class are as defined below:
        INDEX: Stores the index of the node, used for naming and reference purposes
        PROFILE: Stores the profile of the given node. This profile is supplied at the time of initialisation
        STATUS: Status of the node (active/inactive). Every node is node is active during initialisation. A node is 
            inactivated if it has been joined with another node
        TOPHITS: A sequential list of the top hits of a node. The first element in this list is the closest neighbour of
            the node
        CHILDREN: A list of children nodes of a given node. This list is empty during initialisation.
        PARENT: The parent node of a given node.
        AGE: Age of a given node. This variable is used to determine if a node is so old that it top hits need to be re-
            computed
        
    Functions defined in this class are as follows:
        INACTIVATE(): Inactivates the node
        GETAGE(): Returns the age of the node
        SETAGE(): Sets the age of the node to a given value
        GETNODELEVEL(): Returns the level of the node (leaf/internal)
        GETPROFILE(): Returns the profile of the node
        SETPROFILE(): Sets the profile of the node to the given profile
        GETNODESTATUS(): Returns the status of the node (active/inactive)
        GETTOPHITS(): Returns the top hits of the node
        ADDTOPHITS(): Adds a list of top hits to the existing list of top hits and sorts the entire list according to 
                      the neighbour joining score
        CLEARTOPHITS(): Clears the list of top hits
        GETINDEX(): Returns the index of the node
        SETCHILDREN(): Sets a pair of node indices as the children (in the form of a tuple) of the current node
        GETCHILDREN(): Returns the children (in the form of a tuple) of the current node
        SETPARENT(): Sets the given node index as the parent of the current node
        GETPARENT(): Returns the index of the parent node of the current node 
"""
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
