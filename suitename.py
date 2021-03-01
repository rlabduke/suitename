import suiteninit, suitenout
from suitenamedefs import Suite, Residue, Bin, Cluster, Issue, failMessages
from suiteninit import args, bins, MAX_CLUSTERS
from suiteninit import normalWidths, satelliteWidths
from suiteninput import readResidues, readKinemageFile, buildSuites

import sys, os
import numpy as np
from math import cos, pi

#                           suitename.py
# ***************************************************************
# NOTICE: This is free software and the source code is freely
# available. You are free to redistribute or modify under the
# conditions that (1) this notice is not removed or modified
# in any way and (2) any modified versions of the program are
# also available for free.
#               ** Absolutely no Warranty **
# Copyright (C) 2007 David C. Richardson
# ***************************************************************

# 0.2.070524 preserve chi-1 and chi, so could preserve eta, theta
# 0.3.070525 general read dangle record for, e.g.,  eta, theta
# 0.3.070628 triage reports zeta-1, epsilon-1, delta-1,... Ltriage codes
# 0.3.070803 notes: rearranged suitenhead.h/janesviews ...
#                  put something in to say what veiws mean
#  what are masters e and d ????
# 0.3.070919 3g wannabe (tRNA TpseudoUC loop)
# 0.3.110606 range of delta updated by S.J.
# 01/07/2014 S.J. updated so that it can take input with alternate conformations, *nd by default will calculate the suite for altA
# 09/18/2014 S.J. updated so that suitename will ignore DNA residues


version = "suitename.0.6.012521"
dbCounter = 0
dbTarget = 10000  # triggers extra output on this suite

# A collection of variables used for output
class OutNote:
    pass

outNote = OutNote()
outNote.version = version
outNote.comment = ""
outNote.wannabes = 0
outNote.outliers = 0


# ***main()******************************************************************
def main():
    global dbCounter  # for debugging KPB 210222

    # 1. read the input
    if args.infile != "":
        inFile = open(args.infile)
#    elif sys.gettrace() is not None: # how to detect debugger present
    else:
        inFile = sys.stdin

    if args.suitein or args.suitesin:
        suites = readKinemageFile(inFile)
        if len(suites) == 0:
            sys.stderr.write("read no suites: perhaps wrong type of kinemage file\n")
            sys.exit(1)
    else:
        residues = readResidues(inFile)
        if len(residues) == 0:
            sys.stderr.write("read no residues: perhaps wrong alternate code\n")
            sys.exit(1)
        suites = buildSuites(residues)
        suites = suites[:-1]

    # 2. process the suites
    for s in suites:
        if not s.validate():
            if args.test:
                sys.stderr.write(f"! failed validation: {s.pointID}\n")
            suitenout.write1Suite(
                s, bins[13], bins[13].cluster[0], 0, 0, " tangled ", "", "", "", ""
            )
            continue

        # At this point we have a complete suite
        bin, issue, text, pointMaster = evaluateSuite(s)
        pointColor = "white"
        if bin is None:
            s.cluster = bins[0].cluster[0]
            bins[0].cluster[0].count += 1
            suitenout.write1Suite(
                s, bins[0], bins[0].cluster[0], 0, 0, text, issue, "", "", ""
            )
        else:
            memberPack = membership(bin, s)
            (cluster, distance, suiteness, notes, comment,
                pointMaster, pointColor) = memberPack
            suitenout.write1Suite(
                s, bin, cluster, distance, suiteness, notes, issue, comment, 
                pointMaster, pointColor)
            s.suiteness = suiteness
            s.distance = distance
            s.notes = notes
        s.pointMaster = pointMaster
        s.pointColor = pointColor
        dbCounter += 1

    finalStats()
    suitenout.writeFinalOutput(suites, outNote)


# *** evaluateSuite and its tools ***************************************

# Boundaries of various angle ranges
epsilonmin = 155
epsilonmax = 310  # 070130
delta3min = 60
delta3max = 105  #  changed by S.J. on 06/06/2011
delta2min = 125
delta2max = 165
gammapmin = 20
gammapmax = 95  # max 070326
gammatmin = 140
gammatmax = 215  # max 070326
gammammin = 260
gammammax = 335  # max 070326
alphamin = 25
alphamax = 335
betamin = 50
betamax = 290
zetamin = 25
zetamax = 335

# triage table for yes-no angles:
# each of these filters, applied to a suite, will provide a
# true or false answer as to whether this angle is in a reasonable range.
# pointMaster is for grouping points in kinemage display

# data per line: angle index, min, max, code, text, pointMaster
triageFilters = {
    "epsilon": (2, epsilonmin, epsilonmax, Issue.EPSILON_M, "e out", "E"),
    "alpha":   (4, alphamin, alphamax, Issue.ALPHA, " a out", "T"),
    "beta":    (5, betamin, betamax, Issue.BETA, " b out", "T"),
    "zeta":    (3, zetamin, zetamax, Issue.ZETA_M, " z out", "T"),
}

def triage(selector, suite):
    # if angle lies outside the acceptable range, triage immediately
    filter = triageFilters[selector]
    index, min, max, failCode, failText, pointMaster = filter
    if suite.angle[index] < min or suite.angle[index] > max:
        return False, failCode, failText, pointMaster
    else:
        return True, None, None, ""


# The more complex angles are handled by a "sieve".
# A sieve will determine whether an angle is within one of several ranges
# and provide an appropriate code indicating the range.
# This is handled by the sift() function.
sieveDelta = (
    (delta3min, delta3max, 3),
    (delta2min, delta2max, 2),
)

sieveGamma = (
    (gammatmin, gammatmax, "t"),
    (gammapmin, gammapmax, "p"),
    (gammammin, gammammax, "m"),
)

def sift(sieve, angle, failCode):
    for filter in sieve:
        min, max, code = filter
        if min <= angle <= max:
            return code, "", ""
    failMessage = failMessages[failCode]
    return None, failCode, failMessage


def evaluateSuite(suite):
    global bins

    # The order of triage operations, though it may seem arbitrary,
    # was carefully chosen by the scientists.
    ok, failCode, notes, pointMaster = triage("epsilon", suite)
    if not ok:
        return None, failCode, notes, pointMaster

    # Angles with several meaningful ranges:
    # for each angle, find out which range it lies in, or none
    # this becomes a selector to help choose a bin
    puckerdm, failCode, notes = sift(sieveDelta, suite.deltaMinus, Issue.DELTA_M)
    if not puckerdm:
        return None, failCode, notes, "D"

    puckerd, failCode, notes = sift(sieveDelta, suite.delta, Issue.DELTA)
    if not puckerd:
        return None, failCode, notes, "D"

    gammaname, failCode, notes = sift(sieveGamma, suite.gamma, Issue.GAMMA)
    if not gammaname:
        return None, failCode, notes, "T"

    ok, failCode, notes, pointMaster = triage("alpha", suite)
    if not ok:
        return None, failCode, notes, pointMaster

    ok, failCode, notes, pointMaster = triage("beta", suite)
    if not ok:
        return None, failCode, notes, pointMaster

    ok, failCode, notes, pointMaster = triage("zeta", suite)
    if not ok:
        return None, failCode, notes, pointMaster

    # We have pass the test: now use this information to select a bin
    bin = bins[(puckerdm, puckerd, gammaname)]
    # bins is an associated dictionary indexed by the triplet of three angle classifiers
    # each unique triplet of classifiers selects one unique bin, for a total of 12 bins.
    return bin, None, None, ""


# ***membership()***************************************************************

# cluster membership:
# given the bin, we are looking for the correct cluster
def membership(bin, suite):
    matches = np.full(MAX_CLUSTERS, 999.9)
    matchCount = 0
    comment = ""
    pointMaster = ""
    pointColor = "white"

    lDominant = bin.dominant > 0
    if lDominant:  # this bin has a dominant cluster, note it
        dominantJ = bin.dominant
        domCluster = bin.cluster[bin.dominant]

    # find the closest cluster
    # search every cluster in the bin except cluster 0, which is for outliers
    closestD = 999
    for j, c in enumerate(bin.cluster[1:], 1):
        if c.status == "wannabe" and args.nowannabe:
            continue
        distance = hyperEllipsoidDistance(
            suite.angle, bin.cluster[j].angle, 4, normalWidths
        )
        if distance < closestD:
            closestD = distance
            closestJ = j
            closestCluster = c
        matches[j] = distance
        if distance < 1:  # suite could be a member of this cluster
            matchCount += 1

    if matchCount == 1:
        theCluster = closestCluster
        situation = "1-only-one"

    elif matchCount > 1 and not lDominant:
        # dominant cluster is not a possible cluster
        # just output than minimum distance match
        theCluster = closestCluster
        situation = "{matchCount}-None-dom"

    elif matchCount > 1:  # and lDominant
        # find the closest cluster that is not the dominant cluster
        closestNonD = 999
        for j, c in enumerate(bin.cluster[1:], 1):
            if c.status == "wannabe" and args.nowannabe:
                continue
            if matches[j] < closestNonD and c.dominance != "dom":
                closestNonD = matches[j]
                closestJ = j
                theCluster = c

        if theCluster.dominance == "sat":
            # We need to distinguish carefully whether our suite
            # is in the dominant or satellite cluster
            theCluster, closestJ, situation = domSatDistinction(
                suite, domCluster, theCluster, matches, matchCount
            )
        else:
            if matches[dominantJ] < matches[closestJ]:
                closestJ = dominantJ
                theCluster = domCluster
            situation = f"{matchCount}-not-sat"
    else:
        # no match, it's an outlier
        closestJ = 0
        theCluster = closestCluster
        situation = f"outlier distance {closestD:.3}"
        outNote.outliers += 1
        pointMaster = "O"
        pointColor = "white"

    # final computation of suiteness
    # this time we use all 7 dimensions
    if dbCounter >= dbTarget and dbCounter <= dbTarget + 1:  # KPB debug
        print(suite.pointID)
        print(suite.angle)
        print(theCluster.name)
      
    distance = hyperEllipsoidDistance(suite.angle, theCluster.angle, 7, normalWidths)
    # this calculation can assign or deassign a cluster
    if distance <= 1:
        suiteness = (cos(pi * distance) + 1) / 2
        if suiteness < 0.01:
            suiteness = 0.01
    else:
        if closestJ != 0:
            # 7D distance forces this suite to be an outlier
            # so we deassign it here
            closestJ = 0
            comment = f"7D dist {theCluster.name}"
            if theCluster.status == "wannabe":
                comment += " wannabe"
        theCluster = bin.cluster[0]  # outlier
        suiteness = 0

    theCluster.count += 1
    suite.cluster = theCluster
    if theCluster.status == "wannabe" and not args.nowannabe:
        outNote.wannabes = 1  # once set, stays set
    pointColor = 0  # will be handled later!!
    if args.test:
        print(" [suite: %s %s 4Ddist== %f, 7Ddist== %f, suiteness==%f] \n" % \
                (theCluster.name, suite.pointID[:11], closestD, distance, suiteness))
    return theCluster, distance, suiteness, situation, comment, pointMaster, pointColor


def domSatDistinction(suite, domCluster, satCluster, matches, matchCount):
    # if dotproducts both positive, then inbetween
    #      p                  p
    #  dom/___sat  and   dom___\sat
    closestCluster = satCluster
    closestJ = satCluster.ordinal
    dominantJ = domCluster.ordinal

    # use vector properties of numpy.array to determine difference vectors
    domToPoint = domCluster.angle - suite.angle
    satToPoint = satCluster.angle - suite.angle
    domToSat = domCluster.angle - satCluster.angle
    satToDom = -domToSat

    dps = narrowDotProduct(domToPoint, domToSat, 4)
    spd = narrowDotProduct(satToPoint, satToDom, 4)

    if dps > 0 and spd > 0:
        # the trickiest case: point is between dom and sat
        domWidths = normalWidths.copy()
        if args.satellites:
            satWidths = satelliteWidths.copy()
        else:
            satWidths = normalWidths.copy()
        if satCluster.satelliteInfo is not None:
            modifyWidths(domWidths, satWidths, satCluster.satelliteInfo)
        disttodom = hyperEllipsoidDistance(suite.angle, domCluster.angle, 4, domWidths)
        disttosat = hyperEllipsoidDistance(suite.angle, satCluster.angle, 4, satWidths)
        if disttodom < disttosat:
            closestJ = dominantJ
            closestCluster = domCluster
        situation = f"{matchCount}-BETWEEN-dom-sat({disttodom:7.3}|{disttosat:7.3})"
        # else the satellite cluster remains the chosen cluster

    else:
        # the point is not in between
        # just assign by closest standard distance evaluation
        if matches[dominantJ] < matches[closestJ]:
            closestJ = dominantJ
            closestCluster = domCluster
        situation = f"{matchCount}-OUTSIDE-dom-sat"

    return closestCluster, closestJ, situation


# *** Gathering some statistics *********************************************

def finalStats():
    for bin in bins.values():
        for c in bin.cluster:
            bin.count += c.count


# *** The fancy math ********************************************************

# This variable was experimental but we have settled on 3:
power = 3

def hyperEllipsoidDistance(suiteAngles, clusterAngles, nAngles, widthArray):
    global dbCounter
    if nAngles == 4:
        workRange = range(2, 6)
    else:
        workRange = range(1, 8)

    summation = 0
    for k in workRange:
        delta = abs(suiteAngles[k] - clusterAngles[k])
        delta = delta / widthArray[k]
        delToPower = pow(delta, power)
        summation = summation + delToPower
        if dbCounter >= dbTarget and nAngles > 4:  # KPB debug 120221
            sys.stderr.write("db=%3d, k=%d, del=%8.4f, delpower=%10.6f, dpower=%10.6f\n" % 
                        (dbCounter, k, delta, delToPower, summation) )
    result = pow(summation, 1 / power)
    if dbCounter == dbTarget and nAngles > 4:
        sys.stderr.write("final = %7.3f\n" % result)
    return result


# The narrow dot product involves only a subset of the dimensions,
# either 4 or 6. In practice, only 4 is in use.
def narrowDotProduct(a, b, nAngles):
    if nAngles == 4:
        return np.dot(a[2:6], b[2:6])
    else:
        return np.dot(a[1:8], b[1:8])


def modifyWidths(dom, sat, satInfo):
    for m in range(9):
        if satInfo.satelliteWidths[m] > 0:
            sat[m] = satInfo.satelliteWidths[m]
        if satInfo.dominantWidths[m] > 0:
            dom[m] = satInfo.dominantWidths[m]


def showHelpText():
    sys.stderr.write(
        f"""
Version {version}
suitename -flags <stdin >stdout
  or
suitename inputfile -flags >stdout
output flags: [ -report || -string || -kinemage ]
default:  -report

input flags: [ -residuein || -suitein  ]
flags: [ -residuein [ -pointIDfields # ] ] default#=={args.pointidfields}
 OR 
flags: [ -suitein [ -anglefields # ] ]   default#=={args.anglefields}
defaults: -residuein  -pointIDfields {args.pointidfields}

The -residuein format:
label:model:chain:number:ins:type:alpha:beta:gamma:delta:epsilon:zeta
if the file has alternate conformations, then use both -pointIDfields 
    and -altIDfield # to specify the number of pointID fields and which field
    contains the altID
use -altIDval <altID> to specify which alternate conformation to calculate 
    suite for. By default calculated for alt A

-suitein takes a kinemage format,  and uses records from 
    @ballists and/or @dotlists in this format:
{{ptID}} [chi] deltam epsilon zeta alpha beta gamma delta [chi] 
    @dimension in the file, if present, overrides -anglefields

Note dangle trick to make theta,...,eta suites directly

flag: -report [ -chart ]
 suites in order of input, suiteness summary at end
( -chart : NO summary at end, for MolProbity multichart)

flag: -string [-nosequence] [-oneline] 
 3 character per suite string in order of input
    20 per line, ptID of n*20th at end of line
  flag: -nosequence
    only suite names, no Base sequence character
  flag: -oneline
    string all one line, no point IDs

flag: -kinemage
 kinemage of clusters grouped by pucker,pucker ... 
 group {{delta,delta}},subgroup {{gamma}},list {{cluster name}}
  flag: -etatheta or -thetaeta
    kinemage labels theta,eta instead of chi-1,chi

flag: -satellite
  use special general case satellite widths
flag: -nowannabe   
  never assign suites to wannabe clusters
Note: any DNA residues found in the input will be ignored.
""")

main()
