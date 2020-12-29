from suitenamedefs import Suite, Residue, Bin, Cluster
from suiteninit import BuildTheBins
from suiteninput import readResidues, buildSuites
from suitenutil import hyperEllipsoidDistance

import math
import argparse
from enum import Enum

#                           suitename.c                         
#***************************************************************			
# NOTICE: This is free software and the source code is freely   
# available. You are free to redistribute or modify under the   
# conditions that (1) this notice is not removed or modified    
# in any way and (2) any modified versions of the program are   
# also available for free.                                      
#               ** Absolutely no Warranty **                    
# Copyright (C) 2007 David C. Richardson                        
#***************************************************************


#include "suitenscrt.h"
#include "suitenutil.h"
#include "suiteninit.h"
#include "suiteninpt.h"
#include "suitenout.h"

#0.2.070524 preserve chi-1 and chi, so could preserve eta, theta 
#0.3.070525 general read dangle record for, e.g.,  eta, theta 
#0.3.070628 triage reports zeta-1, epsilon-1, delta-1,... Ltriage codes 
#0.3.070803 notes: rearranged suitenhead.h/janesviews ...  
#                  put something in to say what veiws mean 
#  what are masters e and d ???? 
#0.3.070919 3g wannabe (tRNA TpseudoUC loop) 
#0.3.110606 range of delta updated by S.J. 
# 01/07/2014 S.J. updated so that it can take input with alternate conformations, *nd by default will calculate the suite for altA
# 09/18/2014 S.J. updated so that suitename will ignore DNA residues

# reasons why a suite may fail to be classified:
Issue = Enum('Issue', 'DELTA_M EPSILON_M ZETA_M ALPHA BETA GAMMA DELTA')
reasons = {
  Issue.DELTA_M:    " delta-1",
  Issue.EPSILON_M:  " epsilon-1",
  Issue.ZETA_M:     " zeta-1",
  Issue.ALPHA:      " alpha",
  Issue.BETA:       " beta",
  Issue.GAMMA:      " gamma",
  Issue.DELTA:      " delta"
}

def parseCommandLine():
  parser = argparse.ArgumentParser()
  # input styles
  parser.add_argument("residuein")
  parser.add_argument("residuesin")
  parser.add_argument("suitein")
  parser.add_argument("suitesin")
  
  # output styles
  outputStyle = parser.add_mutually_exclusive_group()
  outputStyle.add_argument("report")  # the default format
  outputStyle.add_argument("string")
  outputStyle(add_argument("kinemage"))
  parser.add_argument("chart")     # a modifier to report
  
  # additional options
  parser.add_argument("satellites")
  parser.add_argument("wannabes")
  parser.add_argument("nosequence")

  # numerical options
  parser.add_argument("-v", "--verbosity", type=int,)
  parser.add_argument("ptID", "pointIDfields", type=int)
  parser.add_argument("altIDfield", type=int)
  parser.add_argument("angles", type=int)

  # now actually parse them
  args = parser.parse_args()
  return args


#***main()******************************************************************
def main(argc, argv):
  global args, bins

  args = parseCommandLine()
  bins = buildTheBins()
  # initializations()

  LOK=1
  suiteness=0
  distance=0

  ptcolor[0]  = '\0'  #default is no point color
  ptmaster[0] = '\0'  #default is no point master
  version = "suitename.0.5.122920"

  if args.suitein or args.suitesin:
    # a kinemage format,not yet supported
    pass
  else:
    residues = readResidues()
    suites = buildSuites(residues)


# Boundaries of various angle ranges
epsilonmin = 155; epsilonmax = 310  # 070130
delta3min  =  60; delta3max  = 105  #  changed by S.J. on 06/06/2011 
delta2min  = 125; delta2max  = 165
gammapmin  =  20; gammapmax  =  95  # max 070326
gammatmin  = 140; gammatmax  = 215  # max 070326
gammammin  = 260; gammammax  = 335  # max 070326
alphamin   =  25; alphamax   = 335
betamin    =  50; betamax    = 290
zetamin    =  25; zetamax    = 335

# triage table for yes-no angles:
# each of these filters, applied to a suite, will provide a
# true or false answer as to whether this angle is in a reasonable range.

# data per line: angle index, min, max, code, text
triageFilters =(
  (2, epsilonmin, epsilonmax, Issue.EPSILON_M, "e out"), 
  (4, alphamin, alphamax, Issue.ALPHA, " a out"), 
  (5, betamin, betamax, Issue.BETA, " b out"), 
  (3, zetamin, zetamax, Issue.ZETA_M, " z out"), 
)

# The more complex angles are handled by a "sieve".
# A sieve will determine whether an angle is within one of several ranges
# and provide an appropriate code indicating the range.
sieveDelta = (
# This is handled by the sift() function.
  (delta3min, delta3max, 3),
  (delta2min, delta2max, 2),
)

sieveGamma = (
  (gammatmin, gammatmax, 't'),
  (gammapmin, gammapmax, 'p'),
  (gammammin, gammammax, 'm'),
)


def sift(sieve, angle, failMessage, failCode):
  for filter in sieve:
    min, max, code = filter
    if angle >= min and angle <= max:
      return code, "", ""
  return None, failMessage, failCode


def evaluateSuite(suite):
  pass

  for filter in triageFilters:
    index, min, max, failCode, failText=filter
    if angle[index] < min  or angle[index]>max:
      sour = failText
      triageCode =failCode
      ok = False
      return None, failCode, failText

  puckerdm, failCode, sour = sift(sieveDelta, deltam, Issue.DELTA_M, "'D'",)
  if not puckerdm:
    return None, failCode, sour

  puckerd, failCode, sour = sift(sieveDelta, delta, Issue.DELTA, "'D'",)
  if not puckerd:
    return None, failCode, sour

  gammaname, failCode, sour = sift(sieveGamma, gamma, Issue.GAMMA, "'T'")
  if not gammaname:
    return None, failCode, sour

  # THE BIN IS SELECTED. 
  bin = bins[(puckerdm, puckerd, gammaname)]
  # bins as an associated dictionary indexed by the triplet of three angle classifiers
  # each unique triplet of classifiers selects one unique bin, for a total of 12 bins.
  
  return bin, None, None


def ():
  pass


# cluster membership:
# i is the bin, we are looking for the correct cluster
def membership(suite, i):
  bin = bins[i]
  if bin.dominant > 0:
    lDominant = True

  # find the closest cluster
  j = 1
  for c in bin.clusters:
    if c.status == wannabe  and not args.wannabe:
      j = j + 1; continue
    coordWidths=setUpCoordWidths(bin,j,special=False)
    distance=hyperEllipsoidDistance(suite.angles, bin.cluster[j].angles, 4, coordWidths)
    if distance < closestD:
      closestD = distance
      closestJ = j
    matches[i][j] = distance
    if distance < 1: # suite could be a member of this cluster
      matchCount += 1
    j = j + 1

  # find the next closest cluster
  for c in bin.clusters:
    if c.status == wannabe  and not args.wannabe:
      j = j + 1; continue
    coordWidths=setUpCoordWidths(bin,j,special=False)
    distance=hyperEllipsoidDistance(suite.angles, bin.cluster[j].angles, 4, coordWidths)
    if distance < nextClosestD and j!=closestJ:
      nextClosestD = distance
      nextClosestJ = j
    matches[i][j] = distance
    if distance < 1: # suite could be a member of this cluster
      matchCount += 1
    j = j + 1
  return closestJ, closestD, suiteness, sour   #, issue, comment pointMaster, pointColor



  #    #NOTE: for consistency, edit defaults text in usageout() 
# orphan: pre initializing command line option variables

  # Lreportout = 1 # suite by suite suiteness report, & summary 
  # Lchart = 0     # summary-less report for MolProbity multichart 070521
  # Ldangle = 0  # read straight dangle records  070525
  # Lsourout = 0 # extra info in kinemage ptIDs. optional as of 070524
  # Letatheta = 0 # theta,eta instead of chi-1,chi kinemage angles 070524
  # Lkinemageout=0 # kinemage of the clusters 
  # Lstringout = 0 # 3 char string instead of cluster kinemage 
  # Lsuitesin=0
  # Lresiduesin=1
  # NptIDfields=6  #for dangle residue input
  # altIDfield=0 # by default, no altID field - S.J. 01/07/2014 
  # altID[0]='\0'; altID[1]='\0'
  # Nanglefields=9 #for 9D kinemage edited suite input
  # Lnewfile=0
  # Lhelpout=0
  # LNptIDfields=0 # initializing the four variables to check correct usage of altIDfield - S.J. 01/07/2014
  # LaltIDfield=0
  # Luseincorrect=0 
  # LaltID=0
  # Lchangesout=0
  # NameStr[0] = '\0'
  # Lgeneralsatw = 0 #flag for special general case satellite widths 070328
  # Lwannabe = 1# -wannabe input flag 070429, default 070525, else -nowannabe

  #  Lsequence = 1 #output 1 letter Base code sequence part of string 070409
  # Loverlap = 0  #overlap string lines, e.g. 1-20, 11-30, 21-40, ... 070409
  # Loneline = 0  #string output all as oneline 070409
