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
    
    def __init__(self, name, satelliteWidths, dominantWidths):
        self.name = name
        self.satelliteWidths = satelliteWidths
        self.dominantWidths = dominantWidths


class Residue:
    '''
    # A residue as normally read in, consisting of its six dihedral angles
    Used only briefly as input.
    '''
    pointIDs = []
    base = " "    # A, C, G, U, ...
    angle = np.empty(0)  # will have 6 or 7 elements

    def __init__(self, ID, base, angles):
        self.pointIDs = ID
        self.base = base
        self.angle = angles
        self.unpackAngles()

    def unpackAngles(self):
      pass
    
    # nicknames: for ease of reading the code, each angle is given
    # a meaningful alias. Here they are:
    # 0   alpha
    # 1   beta
    # 2   gamma
    # 3   delta
    # 4   epsilon
    # 5   zeta
    # 7   chi
    @property
    def alpha(self):
        return self.angle[0]

    @alpha.setter
    def alpha(self, value):
        self.angle[0] = value

    @property
    def beta(self):
        return self.angle[1]

    @beta.setter
    def beta(self, value):
        self.angle[1] = value

    @property
    def gamma(self):
        return self.angle[2]

    @gamma.setter
    def gamma(self, value):
        self.angle[2] = value

    @property
    def delta(self):
        return self.angle[3]

    @delta.setter
    def delta(self, value):
        self.angle[3] = value

    @property
    def epsilon(self):
        return self.angle[4]

    @epsilon.setter
    def epsilon(self, value):
        self.angle[4] = value

    @property
    def zeta(self):
        return self.angle[5]

    @zeta.setter
    def zeta(self, value):
        self.angle[5] = value

    @property
    def chi(self):
        return self.angle[6]

    @chi.setter
    def chi(self, value):
        self.angle[6] = value


class Suite:
    '''
    The set of angles forming the linkage BETWEEN residues.
    This is the core data structure used in most operations of the program.
    '''    
    pointID = ()
    base = " "    # A, C, G, U, ...
    angle = np.empty(0) # will become an np.array of 9 angles
    
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
        if angles is None:
          self.angle = np.full(9, 0.0) 
        else:
          self.angle = angles
        self.suiteness = 0.0
        self.distance = 0.0
        self.notes = ""

    def validate(self):
        # make sure that angles deltaMinus through delta are reasonable
        for i in range(1, 8):
            if self.angle[i] < 0 or self.angle[i] > 360:
                return False
        return True

    # nicknames: for ease of reading the code, each angle is given
    # a meaningful alias. Here they are:
    # 0   chiMinus
    # 1   deltaMinus
    # 2   epsilon
    # 3   zeta
    # 4   alpha
    # 5   beta
    # 6   gamma
    # 7   delta
    # 8   chi

    @property
    def chiMinus(self):
        return self.angle[0]

    @chiMinus.setter
    def chiMinus(self, value):
        self.angle[0] = value

    @property
    def deltaMinus(self):
        return self.angle[1]

    @deltaMinus.setter
    def deltaMinus(self, value):
        self.angle[1] = value

    @property
    def epsilon(self):
        return self.angle[2]

    @epsilon.setter
    def epsilon(self, value):
        self.angle[2] = value

    @property
    def zeta(self):
        return self.angle[3]

    @zeta.setter
    def zeta(self, value):
        self.angle[3] = value

    @property
    def alpha(self):
        return self.angle[4]

    @alpha.setter
    def alpha(self, value):
        self.angle[4] = value

    @property
    def beta(self):
        return self.angle[5]

    @beta.setter
    def beta(self, value):
        self.angle[5] = value

    @property
    def gamma(self):
        return self.angle[6]

    @gamma.setter
    def gamma(self, value):
        self.angle[6] = value

    @property
    def delta(self):
        return self.angle[7]

    @delta.setter
    def delta(self, value):
        self.angle[7] = value

    @property
    def chi(self):
        return self.angle[8]

    @chi.setter
    def chi(self, value):
        self.angle[8] = value



