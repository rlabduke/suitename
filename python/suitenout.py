from suitenamedefs import Bin, Cluster
from suitename import Issue, reasons

import sys
from enum import Enum
import numpy as np
from numpy import array

outFile = sys.stdout   # future: might provide a way to change this


def writeString1Suite(suite):
  if args.sequence:
    basestring = suite.base
  else:
    basestring = ":"
  outfile.write(f"{cluster.name}{basestring}")


def writeReport1Suite(suite, bin, cluster, sour, distance, suiteness, pointMaster, pointColor, issue, comment):  # ? LComment, Ltriage
  # 1. write one line of output for this suite
  if issue:
    reason = Reasons[issue]
  elif comment:
    reason = comment
  if cluster.wannabe:
    stray = " wannabe"
  output = f"{suite.point_id} {bin.name} {cluster.name} {suiteness:5.3} {reasons[issue]} {comment}"
  outfile.write(output)

  # 2. gather statistics
  reportCountAll += 1
  
  if bin.ordinal == 0:
    trigCountAll += 1
  elif bin.ordinal < 13:
    suitenessSumAll += suiteness
    binSuiteCountAll += 1

  if cluster.ordinal == 0:
    cluster.suitenessCount[j] += 1
  else:
    cluster.suitenessSum += suiteness
    # report in statistical baskets at intervals of 0.1:
    # everything from 0 to 0.1 goes in bucket 1
    # ... everything from 0.9 to 1.o goes into bucket 10
    if suiteness == 0:
      bucket = 0
    else:
      bucket = 1 + int(suiteness * 10)
    cluster.suitenessCount[bucket] += 1
  pass



