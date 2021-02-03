import suiteninit, suitenout
from suitenamedefs import Suite, Residue, Bin, Cluster, Issue, failMessages
from suiteninit import args, bins, MAX_CLUSTERS
from suiteninit import normalWidths, satelliteWidths
from suiteninput import readResidues, readKinemageSuites, buildSuites
from suitenutil import hyperEllipsoidDistance

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
    if args.infile != "":
        inFile = open(args.infile)
    elif sys.gettrace() is not None:
        #    inFile = open("C:\\Users\\Ken\\Desktop\\Richardson\\suitename\\3gx5_3h-out.dngl")
        inFile = open("C:\\Users\\Ken\\Desktop\\Richardson\\suitename\\4nLf.suitegeom")
    else:
        inFile = sys.stdin

    suiteness = 0
    distance = 0

    if args.suitein or args.suitesin:
        suites = readKinemageSuites(inFile)
        pass
    else:
        residues = readResidues(inFile)
        if residues is None:
            sys.stderr.write("read no residues: perhaps wrong alternate code\n")
            return
        suites = buildSuites(residues)
        suites = suites[:-1]

    for s in suites:
        if not s.validate():
            if args.test:
                sys.stderr.write(f"! failed validation: {s.pointID}\n")
            suitenout.write1Suite(
                s, bins[13], bins[13].cluster[0], 0, 0, " tangled ", "", "", "", ""
            )
            continue
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
            suitenout.write1Suite(                s,                bin,
                cluster,                distance,                suiteness,
                notes,                issue,                comment,
                pointMaster,                pointColor            )
            s.suiteness = suiteness
            s.distance = distance
            s.notes = notes
        s.pointMaster = pointMaster
        s.pointColor = pointColor
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

# data per line: angle index, min, max, code, text
triageFilters = (
    (2, epsilonmin, epsilonmax, Issue.EPSILON_M, "e out", "E"),
    (4, alphamin, alphamax, Issue.ALPHA, " a out", "T"),
    (5, betamin, betamax, Issue.BETA, " b out", "T"),
    (3, zetamin, zetamax, Issue.ZETA_M, " z out", "T"),
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

    # first, the simple angles:
    # if angle lies outside the acceptable range, triage immediately
    for filter in triageFilters:
        index, min, max, failCode, failText, pointMaster = filter
        if suite.angle[index] < min or suite.angle[index] > max:
            notes = failText
            return None, failCode, failText, pointMaster

    # then, the angles with several meaningful ranges:
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

    # use this information to select a bin
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
    #!! print(f"mindistance={closestD:.3f}")

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
        #    theCluster = bin.cluster[0]
        situation = f"outlier distance {closestD:.3}"
        outNote.outliers += 1
        pointMaster = "O"
        pointColor = "white"
    if theCluster.status == "wannabe":
        outNote.wannabes += 1

    # final computation of suiteness
    # this time we use all 7 dimensions
    distance = hyperEllipsoidDistance(suite.angle, theCluster.angle, 7, normalWidths)
    # this calculation can assign or deassign a cluster
    if distance <= 1:
        suiteness = (cos(pi * distance) + 1) / 2
        if suiteness < 0.01:
            suiteness = 0.01
    else:
        if closestJ != 0:
            # 7D distance forces assignment to be an outlier
            # so we deassign it here
            closestJ = 0
            comment = f"7D dist {theCluster.name}"
        theCluster = bin.cluster[0]  # outlier
        suiteness = 0

    theCluster.count += 1
    suite.cluster = theCluster
    pointColor = 0  # will be handled later!!
    if args.test:
        print(" [suite: %s %s 4Ddist== %f, 7Ddist== %f, suiteness==%f] \n" % \
                (theCluster.name, suite.pointID[:11], closestD, distance, suiteness))
    return theCluster, distance, suiteness, situation, comment, pointMaster, pointColor

" [suite: %s %s power== %4.2f, 4Ddist== %f, 7Ddist== %f, suiteness==%f] \n"


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

    if (
        narrowDotProduct(domToPoint, domToSat, 4) > 0
        and narrowDotProduct(satToPoint, satToDom, 4) > 0
    ):
        # the trickiest case: point is between dom and sat
        domWidths = normalWidths.copy()
        satWidths = satelliteWidths.copy()
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
        # bin = bins[key]
        for c in bin.cluster:
            bin.count += c.count


# *** The fancy math ********************************************************

# This variable was experimental but we have settled on 3:
power = 3


def hyperEllipsoidDistance(suiteAngles, clusterAngles, nAngles, widthArray):
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
    result = pow(summation, 1 / power)
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


main()
