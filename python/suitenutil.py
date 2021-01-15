# This file is no longer in use, until further notice


# This variable has been experimental but we have settled on 3
power = 3


def hyperEllipsoidDistance(suiteAngles, clusterAngles, nAngles, widthArray):
    if nAngles == 4: workRange = range(2,6)
    else:            workRange = range(1,8)

    summation = 0
    for k in workRange:
        delta = abs(suiteAngles[k] - clusterAngles[k])
        delta = delta/widthArray[k]
        delToPower = pow(delta,power)
        summation = summation+delToPower
    result=pow(summation, 1/power)
    return result


def hyperEllipsoidDistance2(suiteAngles, clusterAngles, nAngles, widthArray):
    if nAngles == 4: workRange = range(2,6)
    else:            workRange = range(1,8)

    del2Power = [pow(abs(suiteAngles[k] - clusterAngles[k])/widthArray[k], power) for k in workRange]
    summation = sum(del2Power)
    result=pow(summation, 1/power)
    return result
       

    
    
