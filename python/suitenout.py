from suitenamedefs import Bin, Cluster, Issue, reasons
from suiteninit import args

import sys
from enum import Enum
import numpy as np
from numpy import array

outFile = sys.stdout   # future: might provide a way to change this

reportCountAll = 0
trigCountAll = 0
suitenessSumAll = 0
binSuiteCountAll = 0


def write1Suite(suite, bin, cluster, distance, suiteness, sour, issue, comment, pointMaster, pointColor):
  if args.string:
    string1Suite(suite, cluster)
  elif args.kinemage:
    kinemage1Suite(suite, bin, cluster, sour, distance, suiteness, issue, comment, pointMaster, pointColor)
  else:
    report1Suite(suite, bin, cluster, sour, distance, suiteness, issue, comment)


def writeFinalOutput():
  if args.string:
    pass
  elif args.kinemage:
    kinemageFinal()
  else:
    reportFinal()
  

def string1Suite(suite, cluster):
  if args.sequence:
    basestring = suite.base
  else:
    basestring = ":"
  outFile.write(f"{cluster.name}{basestring}")


def report1Suite(suite, bin, cluster, sour, distance, suiteness, issue, comment):  # ? LComment, Ltriage
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
  outIDs = ':'.join(suite.pointId)
  output = f"{outIDs} {bin.binName} {cluster.name} {float(suiteness):5.3f}{reason}{comment}\n"
  outFile.write(output)

  # 2. gather statistics
  reportCountAll += 1
  
  if bin.ordinal == 0:
    trigCountAll += 1
  elif bin.ordinal < 13:
    suitenessSumAll += suiteness
    binSuiteCountAll += 1

  if cluster.ordinal == 0:
    cluster.suitenessCounts[0] += 1
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


def reportFinal():
  if not args.chart:
      pass


def kinemage1Suite(suite, bin, cluster, sour, distance, suiteness, issue, comment, pointMaster, pointColor):
  pass

def writeFinalOutput():
  if args.string:
    pass
  elif args.kinemage:
    pass
  else:
    writeReportSummary()


def writeReportSummary():
  pass
  #!!

#***** kinemage output format *****************************************************

janesviews = '''
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

END
'''

def kinemageHeader():
  # !! put the header here

   outFile.write("@kinemage 1\n");
   outFile.write("@onewidth\n");
   if(args.etatheta): # 070524
     outFile.write("@dimension {theta} {delta-1} {epsilon-1} {zeta-1} {alpha} {beta} {gamma} {delta} {eta}\n")
   else:
     outFile.write("@dimension {chi-1} {delta-1} {epsilon-1} {zeta-1} {alpha} {beta} {gamma} {delta} {chi}\n")
   outFile.write("@dimminmax 0.000 360.000 0.000 360.000 0.000 360.000 0.000 360.000 0.000\
 360.000 0.000 360.000 0.000 360.000 0.000 360.000 0.000 360.000\n");
   if(Lwannabeout):
     outFile.write("@master {wannabees}\n")
   
   # if(LTepsilon) {outFile.write("@pointmaster 'E' {epsilon bad}\n");}
   # else if(LTdelta || LTdeltam) 
   #                       {outFile.write("@pointmaster 'D' {delta bad}\n");}
   # else if(LTzeta || LTalpha || LTbeta || LTgamma)
   #               {outFile.write("@pointmaster 'T' {various bad}\n");}
   # if(Loutlier)  {outFile.write("@pointmaster 'O' {outliers}\n");}

def binGroupOut(deltaMinus, delta):
  from suitename import bins

  for gamma in ('t', 'p', 'm'):
    bin = bins[(deltaMinus, delta, gamma)]
    if bin.active:  # this bin has data points
      outFile.write("@subgroup {{{bin.name}}} recessiveon \n")
      for cluster in bin.clusters:
        if cluster.count > 0:
          extras = ""
          if cluster.status == "wannabe":
            extras = "master= {wannabees}"
          ballList = ("@balllist {{{} {}}} color= {cluster.color} radius= 1 "
                        "nohilite master= {{data}}{}\n").format(
                        bin.name, cluster.name, cluster.color, extras)
          outFile.write(ballList)
          # !!transfer out?
          ringList = ("@ringlist {{{} {}}} color= {} radius= 10 width= 1 "
                        "nobutton master= {{avsigma}} {}\n").format(
              bin.name, cluster.name, cluster.color, extras)
          outFile.write(ringlist)
  pass
