from suitenamedefs import Suite, Residue
from suiteninit import args

import numpy as np
import math

ALTIDFIELD = 4  # where to find codes for alternatives

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
    
  if NAListA.find(baseCode):   base='A'
  elif NAListG.find(baseCode): base='G'
  elif NAListC.find(baseCode): base='C'
  elif NAListU.find(baseCode): base='U'
  elif NAListT.find(baseCode): base='T'
  elif IgnoreDNAList.find(baseCode):
    return None  # we ignore DNA residues
  else:  
    base='Y'
  return base


def stringToFloat(string):
 try:
  n = float(string)
 except ValueError:
  n = 9999.0  # or maybe math.nan?
 return n


def readResidues(inFile):
  lines = inFile.readlines()
  residues = []
  for line in lines:
    if len(line.strip()) == 0 or line[0] == '#':  # blank or comment line
      continue
    fields = line.split(':')
    ids = fields[:args.pointIDfields]
    baseCode = fields[args.pointIDfields-1]
    angleStrings = fields[args.pointIDfields:]
    if ids[ALTIDFIELD][0] != " " and ids[ALTIDFIELD] != args.altID:
      continue  # lines for the wrong alternative conformation are ignored

    base = findBase(baseCode)
    if not base:    # ignore DNA bases
      continue
    angles = np.array([stringToFloat(s) for s in angleStrings])
    for i in range(len(angles)):
        if angles[i] < 0:
            angles[i] += 360.0

    residue = Residue(ids, base, angles)
    residues.append(residue)
  return residues


def buildSuiteBetweenResidues(r1, r2):
  suite = Suite(r2.pointIDs, r2.base)
  suite.chiMinus = r1.chi
  suite.deltaMinus = r1.delta
  suite.epsilon = r1.epsilon
  suite.zeta = r1.zeta
  suite.alpha = r2.alpha
  suite.beta = r2.beta
  suite.gamma = r2.gamma
  suite.delta = r2.delta
  suite.chi = r2.chi

  suite.gatherAngles()
  return suite


def buildSuiteFirst(r2):
  suite = Suite(r2.pointIDs, r2.base)
  suite.alpha = r2.alpha
  suite.beta = r2.beta
  suite.gamma = r2.gamma
  suite.delta = r2.delta
  suite.epsilon = 999
  suite.zeta = 999
  suite.chiMinus = 999
  suite.deltaMinus = 999

  suite.gatherAngles()
  return suite


def buildSuiteLast(r1):
  suite = Suite((),"")
  suite.chiMinus = r1.chi
  suite.deltaMinus = r1.delta
  suite.epsilon = r1.epsilon
  suite.zeta = r1.zeta
  suite.alpha = 999
  suite.beta = 999
  suite.gamma = 999
  suite.delta = 999
  suite.chi = 999

  suite.gatherAngles()
  return suite


def buildSuites(residues): 
  suites = [buildSuiteFirst(residues[0])]
  for i in range(len(residues) - 1):
    suites.append(buildSuiteBetweenResidues(residues[i], residues[i+1]))
  suites.append(buildSuiteLast(residues[-1]))
  return suites


# A kinemage format, not yet supported
def readSuites():
  pass


