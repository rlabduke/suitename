import numpy as np
from numpy import array


class Bin:
    binName = ""
    cluster = ()
    # a tuple of cluster objects
    dominant = -1

    def __init__(self, name, clusters=()):
        # ordinal is just for readability
        self.binName = name
        self.cluster = clusters
        self.dominant = -1

        for i, c in enumerate(clusters):
            if c.dominance == "dominant":
                self.dominant = i
                break



# bins[1] = Bin("33 p")
# bins[1].cluster = (
#     Cluster("!!", "outlier", "white", "out", None, (0,0,0,0,0,0,0,0,0)),
#         ("1a", 1, "certain", "yellowtint", "dominant", one, 
#             (180, 081.495, 212.250, 288.831, 294.967, 173.990, 053.550, 081.035 ,180)))


class Cluster:
    # intrinsic data:
    LOK = False
    name = ""
    status = ""    # certain, wannabe, triaged, outlier, nothing, incomplete
    clusterColor = "" # kinemage color names
    dominance = ""    # dominant, satellite, ordinary, out, tri, inc
    satelliteInfo = None    # present only if this cluster is a satellite
    angles = () 
    # tuple of 9 angles: chi-1 as 0 and chi as 8
    # the standard 7 angles are indices 1-7
    # gathered statistics:
    suitenessSum = 0
    suitenessCounts = None

    def __init__(self, ordinal, name, status, color, dominance, angles):
        self.name = name
        self.LOK = (name != "!!")
        self.status = status
        self.clusterColor = color
        self.dominance = dominance  
        self.angles = angles
        if self.dominance == "sat":
            self.satelliteInfo = satelliteData[name]
        else:
            self.satelliteInfo = None
        self.suitenessCounts = np.zeros(12)
        self.suitenessSum = 0


class SatelliteInfo:
    # numbers used when suite is between satellite and dominant centers
    name = ""
    satelliteWidths = ()
    # vector of 9 angles
    dominantWidths = ()
    # vector of 9 angles
    doma = ""
    # dominant name-for bookkeeping
    
    def __init__(self, name, doma, satelliteWidths, dominantWidths):
        self.name = name
        self.doma = doma
        self.satelliteWidths = satelliteWidths
        self.dominantWidths = dominantWidths


# Residue: a raw residue as normally represented
class Residue:
    pointIDs = []
    base = " "    # A, C, G, U, ...
    alpha = 0
    beta = 0
    gamma = 0
    delta = 0
    epsilon = 0
    zeta = 0

    angle = np.empty(0) 


    def __init__(self, ID, base, angles):
        self.pointIDs = ID
        self.base = base
        self.angle = angles

# Suite: the set of angles forming the linkage BETWEEN residues
class Suite:
    pointId = ""
    base = " "    # A, C, G, U, ...
    chiMinus = 0
    deltaMinus = 0
    epsilon = 0
    zeta = 0
    alpha = 0
    beta = 0
    gamma = 0
    delta = 0
    chi = 0
    angle = np.empty(0) 
    # dual representation: individual angles are named for clarity
    # array is for convenience of computation
    def __init__(self, ID, base, angles=None):
        self.pointId = ID
        self.base = base
        self.angle = angles
        if angles:
            unpackAngles(self)

    def gatherAngles(self):
        self.angle = array(self.chiMinus, self.deltaMinus, self.epsilon, self.zeta, 
            self.alpha, self.beta, self.gamma, self.delta, self.chi)

    def unpackAngles(self):
        self.chiMinus = self.angle[0]
        self.deltaMinus = self.angle[1]
        self.epsilon = self.angle[2]
        self.zeta = self.angle[3]
        self.alpha = self.angle[4]
        self.beta = self.angle[5]
        self.gamma = self.angle[6]
        self.delta = self.angle[7]
        self.k = self.angle[8]



