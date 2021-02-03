import suiteninit

import numpy as np
from numpy import array
from enum import Enum


# reasons why a suite may fail to be classified:
Issue = Enum('Issue', 'DELTA_M EPSILON_M ZETA_M ALPHA BETA GAMMA DELTA')
reasons = {
  Issue.DELTA_M:    "delta-1",
  Issue.EPSILON_M:  "epsilon-1",
  Issue.ZETA_M:     "zeta-1",
  Issue.ALPHA:      "alpha",
  Issue.BETA:       "beta",
  Issue.GAMMA:      "gamma",
  Issue.DELTA:      "delta"
}

failMessages = {
  Issue.DELTA_M:    "bad deltam",
  Issue.GAMMA:      "g out",
  Issue.DELTA:      "bad delta"
}


# primary (coarse grained) classification of suites
class Bin:
    # permanence properties
    name = ""
    ordinal = 0
    cluster = ()
    # a tuple of cluster objects
    dominant = -1
    # statistics gathered during the run
    count = 0
    active = False

    def __init__(self, ordinal, name, clusters=()):
        self.ordinal = ordinal
        self.name = name
        self.cluster = clusters
        self.dominant = -1
        self.active = False

        for i, c in enumerate(clusters):
            if c.dominance == "dom":
                self.dominant = i
                break


# secondary (fine grained) classification of suite
class Cluster:
    # intrinsic data:
    ordinal = 0       # its place in the bin
    name = ""         # the name of the cluster
    status = ""       # certain, wannabe, triaged, outlier, nothing, incomplete
    clusterColor = "" # kinemage color names
    dominance = ""    # dom, sat, ord, out, tri, inc
    satelliteInfo = None    # present only if this cluster is a satellite
  
    # tuple of 9 angles: chi-1 as 0 and chi as 8
    # the standard 7 angles are indices 1-7:
    angle = () 

    # gathered statistics:
    count = 0  # number of data points found in this cluster
    suitenessSum = 0
    suitenessCounts = None

    def __init__(self, ordinal, name, status, color, dominance, angles):
        self.ordinal = ordinal
        self.name = name
        self.LOK = (name != "!!")
        self.status = status
        self.clusterColor = color
        self.dominance = dominance  
        self.angle = array(angles)
        if self.dominance == "sat":
            self.satelliteInfo = suiteninit.getSatelliteInfo(name)
        else:
            self.satelliteInfo = None
        self.suitenessCounts = np.zeros(12)
        self.suitenessSum = 0


class SatelliteInfo:
    # numbers used when suite is between satellite and dominant centers
    name = ""
    satelliteWidths = ()   # vector of 9 angles
    dominantWidths = ()    # vector of 9 angles
    doma = ""    # dominant name-for bookkeeping
    
    def __init__(self, name, doma, satelliteWidths, dominantWidths):
        self.name = name
        self.doma = doma
        self.satelliteWidths = satelliteWidths
        self.dominantWidths = dominantWidths


# Residue: a raw residue as normally represented
class Residue:
    '''
    A residue as normally read in.
    Used only briefly as input.
    '''
    pointIDs = []
    base = " "    # A, C, G, U, ...
    # The 6 angles:
    alpha = 0
    beta = 0
    gamma = 0
    delta = 0
    epsilon = 0
    zeta = 0
    chi = 0
    angle = np.full(7, 0.0)

    def __init__(self, ID, base, angles):
        self.pointIDs = ID
        self.base = base
        self.angle = angles
        self.unpackAngles()

    def unpackAngles(self):
        self.alpha = self.angle[0]
        self.beta = self.angle[1]
        self.gamma = self.angle[2]
        self.delta = self.angle[3]
        self.epsilon = self.angle[4]
        self.zeta = self.angle[5]
        # for the future:
        # we can accept chi as a seventh angle if provided
        if len(self.angle) > 6:
            self.chi = self.angle[6]
        else:
            self.chi = -431602080.00  #180
            # A preposterous compromise with the past for now
            #self.chi = -180
    

# Suite: the 
class Suite:
    '''
    The set of angles forming the linkage BETWEEN residues.
    This is the core data structure used in most operations of the program.
    '''    
    pointID = ()
    base = " "    # A, C, G, U, ...
    #chiMinus = 0
    deltaMinus = 0
    epsilon = 0
    zeta = 0
    alpha = 0
    beta = 0
    gamma = 0
    delta = 0
    chi = 0
    # dual representation: individual angles are named for clarity
    # array is for convenience of computation
    angle = np.full(9, 0.0) 

    @property
    def chiMinus(self):
        return self.angle[0]

    @chiMinus.setter
    def chiMinus(self, value):
        self.angle[0] = value

    # fields computed during analysis:
    cluster = None  # The cluster to which it is assigned
    suiteness = 0.0
    distance = 0.0
    notes = ""
    pointMaster = ""
    pointColor = ""

    def __init__(self, ID, base, angles=None):
        self.pointID = ID
        self.base = base
        if angles is not None:
          self.angle = angles
          self.unpackAngles()
        self.cluster = None
        self.suiteness = 0.0
        self.distance = 0.0
        self.notes = ""

    def validate(self):
        # make sure that angles deltaMinus through delta are reasonable
        for i in range(1, 8):
            if self.angle[i] < 0 or self.angle[i] > 360:
                return False
        return True

    def gatherAngles(self):
        self.angle = array((self.chiMinus, self.deltaMinus, self.epsilon, self.zeta, 
            self.alpha, self.beta, self.gamma, self.delta, self.chi))

    def unpackAngles(self):
        self.chiMinus = self.angle[0]
        self.deltaMinus = self.angle[1]
        self.epsilon = self.angle[2]
        self.zeta = self.angle[3]
        self.alpha = self.angle[4]
        self.beta = self.angle[5]
        self.gamma = self.angle[6]
        self.delta = self.angle[7]
        self.chi = self.angle[8]



