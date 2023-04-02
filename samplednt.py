#samplednt.py implements survivorship in neutral theory with incomplete sampling
#corresponds to equations 8-15 in the methods
import numpy as np, random, matplotlib.pyplot as plt, math, mpmath, time, warnings,sys
sysmin = sys.float_info.min #used when calculations involve very small numbers

#probability of observing each in a list of [delta t]s
#also takes hazard rates,detection rates, and endpoint - the extent to calculate to
#if normalize is true, then condition probabilities on observing duration>0
#when D gets really small, use mpmath
def probdelt(delt,d,h,endpoint,normalize=False):
    S = list(np.cumprod([1] + [1-i for i in h][:-1])) #survivorship to time t
    f = list(np.multiply(h,S)) #chance of extinction between time t and t+1
    D = list(np.cumprod([1-i for i in d])) #probability of no detections between times 0 and t
    if D[-1] < sysmin:
        print("Calculating slowly with mpmath")
        D = list(np.cumprod([mpmath.mpf(1-i) for i in d])) #allows very small numbers
        D = [float(i) if i > sysmin else i for i in D] #only smallest are mpf
    else:
        print("Calculating quickly without mpmath")
    x = [f[tprime]*D[tprime] for tprime in range(endpoint)][::-1]
    middleterm = np.cumsum(x)[::-1]
    pnotobserved = sum([f[t]*D[t] for t in range(endpoint)])
    pdelt0 = sum([d[t]/(1-d[t])*middleterm[t] for t in range(endpoint)])
    out = []
    for dt in delt:
        if np.isnan(dt): #not detected
            out.append(pnotobserved)
        elif dt == 0:
            out.append(pdelt0)
        else:
            p = sum([d[t+dt]/D[t+dt]*middleterm[t+dt]*(D[t-1] if t > 0 else 1)*d[t] for t in range(endpoint-dt)])
            out.append(p)
    if normalize:
        return([i/(1-pnotobserved-pdelt0) for i in out])
    else:
        return(out)

#probcumulate() takes a vector of probability masses whose length is divisible by divisor
def probcumulate(vec,divisor):
    return([sum(vec[i:(i+divisor)]) for i in range(0,len(vec),divisor)])