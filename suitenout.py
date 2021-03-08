from suitenamedefs import Bin, Cluster, Issue, reasons
from suiteninit import args, bins

import sys
from enum import Enum
import numpy as np
from numpy import array

outFile = sys.stdout  # future: might provide a way to change this

reportCountAll = 0
trigCountAll = 0
suitenessSumAll = 0
binSuiteCountAll = 0


def write1Suite(suite, bin, cluster, distance, suiteness, notes, 
                issue, comment, pointMaster, pointColor):
    if args.string:
        string1Suite(suite, cluster)
    elif args.kinemage:
        kinemage1Suite(suite, bin, cluster, notes, distance, suiteness, 
                       issue, comment, pointMaster, pointColor,
        )
    else:
        report1Suite(suite, bin, cluster, notes, distance, suiteness, issue, comment)


def writeFinalOutput(suites, outNote):
    if args.satellites:
        outNote.comment = " using special general case satellite widths"
    if args.string:
        pass
    elif args.kinemage:
        kinemageFinal(suites, outNote)
    else:
        reportFinal(outNote)


def string1Suite(suite, cluster):
    if args.nosequence:
        basestring = ":"
    else:
        basestring = suite.base
    outFile.write(f"{cluster.name}{basestring}")


def report1Suite(
    suite, bin, cluster, notes, distance, suiteness, issue, comment
):  # ? LComment, Ltriage
    global reportCountAll, trigCountAll, suitenessSumAll, binSuiteCountAll

    # 1. write one line of output for this suite
    reason = ""
    if issue:
        reason = " " + reasons[issue]
    elif comment:
        reason = " " + comment
        comment = ""
    if cluster.status == "wannabe":
        comment = " wannabe"
    outIDs = ":".join(suite.pointID)
    output = (
        f"{outIDs} {bin.name} {cluster.name} {float(suiteness):5.3f}{reason}{comment}\n"
    )
    outFile.write(output)

    # 2. gather statistics
    reportCountAll += 1

    if bin.ordinal == 0:
        trigCountAll += 1
    elif bin.ordinal < 13:
        suitenessSumAll += suiteness
        binSuiteCountAll += 1

    if cluster.ordinal == 0:
        cluster.suitenessCounts[11] += 1
    else:
        cluster.suitenessSum += suiteness
        # report in statistical baskets at intervals of 0.1:
        # everything from 0 to 0.1 goes in bucket 1
        # ... everything from 0.9 to 1.o goes into bucket 10
        if suiteness == 0:
            bucket = 0
        else:
            bucket = 1 + int(suiteness * 10)
        cluster.suitenessCounts[bucket] += 1


def reportFinal(outNote):
    if not args.chart:
        outFile.write(outNote.comment + "\n")
        suitenessAverage(0)
        if bins[1].cluster[1].count > 0:  # Aform 1a    070325
            suitenessAverage(1)
            suitenessAverage(2)


def suitenessAverage(mode):
    """
    Gather statistics on suiteness.
    12 buckets:
      one for suiteness=0
      one each for divisions by tenths from 0 to 1
      one for outliers
    """
    bucket = np.zeros(12, dtype=int)
    sum = 0
    average = 0
    allCount = 0
    excludedCluster = None

    if mode == 1:
        # cluster 1a all by itself
        comment = " A form (1a)"
        cluster = bins[1].cluster[1]
        sum = cluster.suitenessSum
        for k in range(12):
            bucket[k] = cluster.suitenessCounts[k]
    else:
        if mode == 0:
            # all complete suites
            startWith = 0
            comment = "For all"
            outFile.write(
                "Found {} complete suites derived from {} entries\n".format(
                    binSuiteCountAll + trigCountAll, reportCountAll
                )
            )
            outFile.write(
                "{} suites were triaged, leaving {} assigned to bins\n".format(
                    trigCountAll, binSuiteCountAll
                )
            )
        elif mode == 2:
            # all complete suites except cluster 1a
            startWith = 1
            comment = " non-1a  has"
            excludedCluster = bins[1].cluster[1]
        for i in range(1, 13):
            bin = bins[i]
            # the final bin 13 is for the pseudo-suites with incomplete
            # angles and is ignored
            for cluster in bin.cluster[startWith:]:
                # cluster 0 in every bin is for outliers and is ignored
                if cluster is excludedCluster:
                    continue  # mode 2: ignore cluster 1a
                sum += cluster.suitenessSum
                # print(
                #     f"cluster:{cluster.name} {cluster.count:2} cluster sum: {cluster.suitenessSum} sum: {sum}\n"
                # )

                for k in range(12):
                    bucket[k] += cluster.suitenessCounts[k]

    allCount = np.sum(bucket)
    if allCount > 1:
        average = sum / allCount
    else:
        average = 0

    outFile.write(
        "{} {} suites: average suiteness == {:5.3f}\n".format(
            comment, allCount, average
        )
    )
    if mode == 0:
        outFile.write(f"{bucket[11]:6d} suites are  outliers\n")
    outFile.write(f"{bucket[0]:6d} suites have suiteness == 0    \n")
    outFile.write(f"{bucket[1]:6d} suites have suiteness >  0 <.1\n")
    for k in range(2, 10):
        outFile.write(
            "{:6d} suites have suiteness >=.{} <.{}\n".format(bucket[k], k - 1, k)
        )
    outFile.write(f"{bucket[10]:6d} suites have suiteness >=.9    \n")



# ***** kinemage output format *****************************************************

def kinemage1Suite(suite, bin, cluster, notes, distance, suiteness, issue, comment, 
                    pointMaster, pointColor):
    suite.pointMaster = pointMaster
    suite.pointColor = pointColor
    suite.notes = notes


# static text: viewing parameters
janesviews = """
@viewid {d e z}
@zoom 1.00
@zslab 200
@ztran 0
@center 197.500 172.300 178.300
@axischoice 2 3 4
@matrix
0.07196 0.11701 -0.99052 -0.00336 0.99312 0.11707 0.99740 -0.00509 0.07186
@2viewid {zag front}
@2zoom 1.00
@2zslab 200
@2ztran 0
@2center 174.091 194.887 207.768
@2axischoice 4 5 7
@2matrix
0.99508 -0.00018 -0.09905 -0.00135 -0.99993 -0.01172 -0.09904 0.0118 -0.99501
@3viewid {a b g}
@3zoom 1.00
@3zslab 200
@3ztran 0
@3center 175.700 189.600 64.100
@3axischoice 5 6 7
@3matrix
0.99955 0.000101 0.030002 0.0002 0.99995 -0.010012 -0.030001 0.010013 0.9995

"""

# static text: static items in the display
kinemageFrame = """
@group {frame} dominant 
@vectorlist {frame} color= white 
P   0.000   0.000   0.000   5.000   0.000   0.000 
P  35.000   0.000   0.000  40.000   0.000   0.000 
P  80.000   0.000   0.000 160.000   0.000   0.000 
P   0.000   0.000   0.000   0.000   5.000   0.000 
P   0.000  35.000   0.000   0.000  40.000   0.000 
P   0.000  80.000   0.000   0.000 160.000   0.000 
P   0.000   0.000   0.000   0.000   0.000   5.000 
P   0.000   0.000  35.000   0.000   0.000  40.000 
P   0.000   0.000  80.000   0.000   0.000 160.000 
P 200.000   0.000   0.000 280.000   0.000   0.000 
P 320.000   0.000   0.000 360.000   0.000   0.000 
P   0.000 200.000   0.000   0.000 280.000   0.000 
P   0.000 320.000   0.000   0.000 360.000   0.000 
P   0.000   0.000 200.000   0.000   0.000 280.000 
P   0.000   0.000 320.000   0.000   0.000 360.000 
@labellist {XYZ} color= white 
{X}  20.000  -5.000  -5.000 
{X} 380.000  -5.000  -5.000 
{Y}  -5.000  20.000  -5.000 
{Y}  -5.000 380.000  -5.000 
{Z}  -5.000  -5.000  20.000 
{Z}  -5.000  -5.000 380.000 
@labellist {mtp} color= green 
{p}  60.000   0.000   0.000 
{t} 180.000   0.000   0.000 
{m} 300.000   0.000   0.000 
{p}   0.000  60.000   0.000 
{t}   0.000 180.000   0.000 
{m}   0.000 300.000   0.000 
{p}   0.000   0.000  60.000 
{t}   0.000   0.000 180.000 
{m}   0.000   0.000 300.000 

"""

def kinemageFinal(suites, outNote):
    """
    Output the content of a kinemage file
    The 3, 2 order may seem odd, but it is a standard
    """
    kinemageHeader(outNote)
    for deltaMinus in (3, 2): 
      for delta in (3, 2):
          binGroupOut(deltaMinus, delta, suites)
    triaged = bins[0]
    if triaged.count > 0:
        outFile.write("@group {triaged} dominant dimension=9 wrap=360 select off\n")
        binOut(triaged, suites)


def kinemageHeader(outNote):
    """ The invariant portion of a kinemage file """
    outFile.write(f"@text\n {outNote.version}\n {outNote.comment}\n")
    outFile.write("@kinemage 1\n")
    outFile.write("@onewidth\n")
    if args.etatheta or args.thetaeta:  # 070524
        outFile.write(
            "@dimension {theta} {delta-1} {epsilon-1} {zeta-1} {alpha} "
            "{beta} {gamma} {delta} {eta}\n"
        )
    else:
        outFile.write(
            "@dimension {chi-1} {delta-1} {epsilon-1} {zeta-1} {alpha} "
            "{beta} {gamma} {delta} {chi}\n"
        )
    outFile.write(
        "@dimminmax 0.000 360.000 0.000 360.000 0.000 360.000 0.000 "
        "360.000 0.000 360.000 0.000 360.000 0.000 360.000 0.000 360.000 "
        "0.000 360.000\n"
        )
    if outNote.outliers:
      outFile.write("@pointmaster 'O' {outliers}\n")
    if outNote.wannabes:
      outFile.write("@master {wannabees}\n")
    outFile.write(janesviews)
    outFile.write(kinemageFrame)


def formatAngles(angles, separator):
    strings = ["{:7.3f}".format(a) for a in angles]
    out = separator.join(strings)
    return out


def binGroupOut(deltaMinus, delta, suites):
    # If any bin in the group has data, generate a group header
    groupCount = 0
    for gamma in ("p", "t", "m"):
      bin = bins[(deltaMinus, delta, gamma)]
      groupCount += bin.count
    if groupCount > 0:        
      outFile.write(
        f"@group {{{deltaMinus}{delta}}} recessiveon dimension=9"
        " wrap=360 select animate off\n")

    # generate the data
    for gamma in ("p", "t", "m"):
        bin = bins[(deltaMinus, delta, gamma)]
        if bin.count > 0:
          binOut(bin, suites)


def binOut(bin, suites):
    if any([c.count > 0 for c in bin.cluster]):
        outFile.write(f"@subgroup {{{bin.name}}} recessiveon \n")
    for cluster in bin.cluster[1:]:
        # the first cluster, for outliers, will be handled later
        if cluster.count > 0:
            extras = ""
            if cluster.status == "wannabe":
                extras = " master= {wannabees}"
            # display a ball for each for each suite in this cluster
            ballList = (
                "@balllist {{{} {}}} color= {} radius= 1 "
                "nohilite master= {{data}}{}\n"
            ).format(bin.name, cluster.name, cluster.clusterColor, extras)
            # display a ring surrounding the center of the cluster
            ringList = (
                "@ringlist {{{} {}}} color= {} radius= 10 width= 1 "
                "nobutton master= {{avsigma}}{}\n"
            ).format(bin.name, cluster.name, cluster.clusterColor, extras)
            angleList = formatAngles(cluster.angle[1:-1], ' ')
            ringList2 = "{{{} {}}} 180 {} 180\n".format(
                bin.name, cluster.name, angleList
            )
            labelList = (
                "@labellist {{{} {}}} color= {} nobutton "
                "master= {{labels}}{}\n"
            ).format(bin.name, cluster.name, cluster.clusterColor, extras)

            outFile.write(ballList)
            outPoints(bin, cluster, suites, "")
            outFile.write(ringList)
            outFile.write(ringList2)
            outFile.write(labelList)
            outFile.write(ringList2)
            
    # handle outliers if there are any:
    if bin.cluster[0].count > 0:
        cluster = bin.cluster[0]
        ballList = (
            "@balllist {{{} {}}} color= {} radius= 1 "
            "nohilite master= {{data}}\n"
        ).format(bin.name, cluster.name, cluster.clusterColor)
        outFile.write(ballList)
        outPoints(bin, cluster, suites, "")


def outPoints(bin, cluster, suites, extra1):
    for s in suites:
        if s.cluster is cluster:
            if bin.name == "trig":  # the triage bin is especially handled
                extra = f"'{s.pointMaster}'"
            else:
                extra = extra1
            ids = ':'.join(s.pointID)
            line = \
              (f"{{{bin.name} {cluster.name} :D=={s.distance:5.3f}"
               f":S=={s.suiteness:5.3f}: {ids}}} {extra} ,") \
               + formatAngles(s.angle, ",") + "\n"
            outFile.write(line)
   
