from suitenamedefs import Suite, Residue
from suitename import args
import math

# The great variety of codes that may represent each base in the input file
 NAListA = ":ADE:  A:A  : Ar:ATP:ADP:AMP:T6A:1MA:RIA:  I:I  :"
 NAListG = ":GUA:  G:G  : Gr:GTP:GDP:GMP:GSP:1MG:2MG:M2G:OMG: YG: 7MG:YG :"
 NAListC = ":CYT:  C:C  : Cr:CTP:CDP:CMP:5MC:OMC:"
 NAListU = ":URA:URI:  U: Ur:U  :UTP:UDP:UMP:5MU:H2U:PSU:4SU:"
 NAListT = ":THY:  T:T  : Tr:TTP:TDP:TMP:"
 IgnoreDNAList = ": DA: DG: DC: DT:"


# out of the noise, determine the base
def findBase(baseCode):
 if len(baseCode) != 3:
  return 'Z'

 if NaListA.find(baseCode):   base='A'
 elif NaListG.find(baseCode): base='G'
 elif NaListC.find(baseCode): base='C'
 elif NaListU.find(baseCode): base='U'
 elif NaListT.find(baseCode): base='T'
 elif IgnoreDNAList.find(baseCode):
  return None  # we ignore DNA residues
 else:  base='Y'
  return base


def string_to_float(string):
 try:
  n = float(string)
 except ValueError:
  n = 9999.0  # or maybe math.nan?
 return n


def readResidues(inFile):
    lines = inFile.readlines()
    residues = []
    for line in lines:
        fields = line.split(':')
        ids = fields[:args.pointIdFields]
        baseCode = fields[args.pointIdFields]
        angleStrings = fields[args.pointIdFields+1:]

        base = findBase(baseCode)
        if not base:    # ignore DNA bases
         continue
        angles = np.array([stringToFloat(s) for s in angleStrings])

        residue = Residue(ids, base, angles)
        residues.append(residue)
        return residues


def buildSuiteBetweenResidues(r1, r2):
  suite = Suite(r2.PointIDs, r2.base)
  suite.chiMinus = r1.chi
  suite.deltaMinus = r1.delta
  suite.epsilon = r1.epsilon
  suite.zeta = r1.zeta
  suite.alpha = r2.alpha
  suite.beta = r2.beta
  suite.gamma = r2.gamma
  suite.delta = r2.delta

  suite.gather_angles()
  return suite


def buildSuites(residues): 
  suites = []
  for i in range(len(residues) - 1):
    suites.append(buildSuiteBetweenResidues(suites[i], suites[i+1]))
  return suites


# A kin image format, not yet supported
def readSuites():
    pass


