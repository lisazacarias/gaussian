#!/usr/local/lcls/package/python/current/bin/python
################################################################################
# Modified version of John Sheppard's read_xcor_data script
# Reads in a matlab file of xcor data
# Data must have column names posList, ampList, and ampstdList
################################################################################

from __future__ import division

from pylab import *
from scipy import special
import scipy.io as sio
import time
from scipy.optimize import curve_fit
import operator
from numpy import inf
from collections import defaultdict
from sys import argv, exit

NUM_BUCKS = 10

def extract(axdata, column):
    matLst = axdata['data'][column][0][0]
    out = []
    for lst in matLst:
        out.append(lst[0])
    return array(out)

def getBucket(val, step):
    bucket = int(floor(val/step))
    return bucket if bucket < 10 else 9

# There has to be a more efficient way of doing this, I'm just tired
def findMax(data, run):
    max = 0
    maxIdx = run[0]
    for idx in run:
        val = data[idx]
        if val > max:
            max = val
            maxIdx = idx
    return maxIdx

# This can probably be done in preprocessing
# The ever-present struggle between speed and space
def findRun(idx, runs):
    for i, run in enumerate(runs):
        if idx in run:
            return i

def findWidths(peakIdxs, runs):
    widths = []
    for peakIdx in peakIdxs:
        runIdx = findRun(peakIdx, runs)
        if runIdx == 0:
            widths.append((runs[1][0] - peakIdx)*2)
        elif runIdx == len(runs)-1:
            widths.append((peakIdx - runs[-2][-1])*2)
        else:
            widths.append(runs[runIdx+1][0] - runs[runIdx-1][-1])
    return widths

# Modified from StackOverflow
def func(x, *params):
    m = params[0]
    b = params[1]

    y = [m*i + b for i in x]
    
    for i in range(2, len(params), 3):
        ctr = params[i]
        amp = params[i+1]
        wid = params[i+2]
        y = y + gaussian(x, ctr, wid, amp)
    return y

def gaussian(x, mu, sig, amp):
    return amp*np.exp(-np.power(x - mu, 2.) / (2 * np.power(sig, 2.)))
    
def getSlope(x1, y1, x2, y2):
    return (y2-y1)/(x2-x1)
    
# Idea to add a line instead of a really short, fat gaussian was all Ahemd.
# Thanks, yo. You're great.
def findLine(zeros, runs, data, nonZeroRuns):
    x1, y1, x2, y2, m, b = (0, 0, 0, 0, 0, 0)

    # TODO: this seriously needs testing
    # This condition should only be possible if there are peaks on one or both 
    # extremes, or if there is no peak
    if len(zeros) == 1:
        zeroRun = runs[zeros[0]]
        # This should pull out the median index value of the run
        x1 = zeroRun[np.argsort(data[zeroRun])[len(zeroRun)/2]]
        y1 = data[x1]
        return [m, y1]
        
    # 0 shouldn't be possible given that the data is normalized to the lowest 
    # point, so it should be 2+
    # Currently just fitting the first point of the first zero run and the last
    # point of the last zero run. Could make it smarter by adding a sum of step
    # functions, but that seems like overkill
    else:
        zero1 = runs[zeros[0]]
        zero2 = runs[zeros[-1]]
        
        x1 = zero1[np.argsort(data[zero1])[len(zero1)/2]]
        y1 = data[x1]
        
        x2 = zero2[np.argsort(data[zero2])[len(zero2)/2]]
        y2 = data[x2]
        
        m = getSlope(x1, y1, x2, y2)
        
        return [m, y1-m*x1]        

def getPeaks(data, numPeaks, runs):
    # Should be doable in preprocessing
    lenRuns = [len(run) for run in runs]
    
    # User-proofing. Could probably limit input
    numPeaks = numPeaks if numPeaks <= len(runs) else len(runs)
    
    # Would be using linear argpartsort if we were running not 2013 builds.
    # Can you tell I'm bitter?
    ind = np.argsort(array(lenRuns))[-numPeaks:]

    # This is inelegant
    peakIdx = []
    for run in array(runs)[ind]:
        peakIdx.append(findMax(data, run))

    peaks = []
    for idx in peakIdx:
        peaks.append(data[idx])

    max_index, max_value = max(enumerate(data), key=operator.itemgetter(1))

    # Maybe unnecessary precaution to make sure that the max point is used in
    # the fit (a run wouldn't be found if the peak were sufficiently narrow)
    if max_value not in peaks:
        min_index, _ = min(enumerate(peaks), key=operator.itemgetter(1))
        peaks[min_index] = max_value
        peakIdx[min_index] = max_index

    return[peaks, peakIdx]


# Checking for inflection points doesn't work because some data points don't 
# follow the trend line. This groups consecutive data points by bucket.
# Buckets need to be recalculated following an adjustment.
def getRuns(data, step):
    zeros = []
    nonZeroRuns = []
    runs = []
    currRun = []
    currBuck = getBucket(data[0], step)
    
    for idx, point in enumerate(data):
        newBuck = getBucket(point, step)
        if newBuck == currBuck:
            currRun.append(idx)
        else:
            # Plotting the end of a run
            # plt.axvline(x=idx)
            
            # Three points make a curve!
            if len(currRun) > 2:
                runs.append(currRun)
                if currBuck == 0:
                    zeros.append(len(runs)-1)
                else:
                    nonZeroRuns.append(currRun)

            currRun = [idx]
            currBuck = newBuck
            
    # Effectively flushing the cache
    if len(currRun) > 2:
        runs.append(currRun)
        if currBuck == 0:
            zeros.append(len(runs)-1)
            
    return [runs, zeros, nonZeroRuns]
    
# A whole rigmarole to collapse multiple pedestals.
# It assumes that the pedestal is the bucket with the most elements
def adjustData(data, step, normalizedAdjustment):
    bucketCount = defaultdict(int)
    bucketContents = defaultdict(list)
    buckets = [0 for i in xrange(0,len(data))]
    
    for idx,element in enumerate(data):
        bucket = getBucket(element, step)
        bucketCount[bucket] += 1
        bucketContents[bucket] += [idx]
        buckets[idx] = bucket
    
    zeroBucket = max(bucketCount.iteritems(), key=operator.itemgetter(1))[0]
    
    needsAdjustment = False
    
    for idx, bucket in enumerate(buckets):
        if bucket < zeroBucket:
            # Inefficient to set this every time, but eh
            needsAdjustment = True
            # Sets them arbitrarily to the value of the first element in the
            # zero bucket, to eliminate the double pedestal
            data[idx] = data[bucketContents[zeroBucket][0]]
    
    if needsAdjustment:
        normalizedAdjustment = min(data[bucketContents[zeroBucket]])
        data = data - normalizedAdjustment
        step = max(data) / NUM_BUCKS
    
    return [data, step, normalizedAdjustment]
    
def getGuess(data, step, useZeros):
    # I feel like I should rename this function to something less... runny
    runs, zeros, nonZeroRuns = getRuns(data, step)
                
    peaks, peakIdx = (getPeaks(data, numPeaks, nonZeroRuns) 
                        if not useZeros 
                        else getPeaks(data, numPeaks, runs))
                        
    # This plots my guesses for the peaks
    #for idx in peakIdx:
    #    plt.axvline(x=idx)

    widths = findWidths(peakIdx, runs)
    
    guess = findLine(zeros, runs, data, nonZeroRuns)
    # This plots my guess for the line
    #plt.plot([guess[0]*j + guess[1] for j in xrange(0, len(data))], '--')

    for idx, amp in enumerate(peaks):
        guess += [peakIdx[idx], amp, widths[idx]/4]
        # Plot my initial guesses for the gaussian(s)
        #plt.plot([gaussian(i, peakIdx[idx], widths[idx]/4, amp)
                  #for i in xrange(0,len(data))], '--')
    
    return guess

def plotFit(data, numPeaks, useZeros):

    #Plot the data
    plt.plot(data, '.', marker='o')

    firstAdjustment = min(data)
    normalizedAdjustment = 0

    # Removing the pedestal
    data = data - firstAdjustment
    
    # Define the step size by the number of vertical buckets
    step = max(data) / NUM_BUCKS
    
    data, step, normalizedAdjustment = adjustData(data, step, 
                                                  normalizedAdjustment)
                                                  
    # This prints my vertical buckets
    #for i in xrange(1,NUM_BUCKS):
        #plt.plot([i*step for _ in xrange(0, len(data))])
        
    totalAdjustment = firstAdjustment + normalizedAdjustment
        
    print "adjustment: " + str(totalAdjustment)

    # Note that the format of the guess is a list of the form:
    # [m, b, center_0, amplitude_0, width_0,..., center_k, amplitude_k, width_k]
    # where m and b correspond to the line parameters in y = m*x + b
    # and every following group of three corresponds to a gaussian 
    guess = getGuess(data, step, useZeros)

    x = range(0,len(data))

    popt, pcov = curve_fit(func, x, data, p0=guess,
                           # Someday this feature will be available...
                           # ...When we're no longer running builds from 2013 :P
                           #bounds=(0, [len(data), max(data),len(data)]),
                           maxfev=20000)
                           

    # Print and plot the optmized line fit.
    # Note that popt has the same format as the guess, meaning that the first
    # two parameters are the m and b of the line, respectively
    print "line: " + "m = " + str(popt[0]) + ", b = " + str(popt[1])
    plt.plot([popt[0]*j + popt[1] + totalAdjustment 
             for j in xrange(0, len(data))], '--')
    
    # Print and plot the optimized gaussian fit(s)
    # Again, the first two elements were the line, and each gaussian is a
    # subsequent group of 3 elements (hence starting at index 2 and incrementing
    # by 3)
    for i in xrange(2, len(popt), 3):
        print ("gaussian " + str(i//3) + ": center = " + str(popt[i])
               + ", amplitude = " + str(popt[i+1]) + ", width = "
               + str(popt[i+2]))
        plt.plot([gaussian(j, popt[i], popt[i + 2], popt[i + 1])
                 + totalAdjustment for j in xrange(0, len(data))], '--')

    fit = func(x, *popt)

    plt.plot(fit + totalAdjustment, linewidth=2)

    show()

if __name__ == "__main__":
    try:
        filepath = argv[1]
    except IndexError:
        print "Usage: " + argv[0] + " [path to XCor matlab file]"
        exit()
        
    axdata = sio.loadmat(filepath)

    ampList = extract(axdata, 'ampList')
    
    # TODO: Ask Axel if he needs these
    #posList = extract(axdata, 'posList')
    #ampstdList = extract(axdata, 'ampstdList')

    numPeaks = input("Number of gaussians to fit: ")
    plotFit(ampList, numPeaks, False)
