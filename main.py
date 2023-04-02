#main.py is used for all the main analyses except fitting NT models, which is done in cumulated.py
from scipy import stats
import time, numpy as np, matplotlib.pyplot as plt,math, pandas as pd
from markov import transmat, hazardandomegastar
from samplednt import probdelt

#surv() takes a list of empirical geological durations and outputs the number of species surviving at least X timesteps up to the max duration
def surv(durations):
    return([sum([v>i for v in durations]) for i in range(int(max(durations)))])

#main loop
if __name__ == '__main__':
    np.random.seed(1)

    #parameters
    #unitlength = 10000 #a graptolite unit has this expected lifespan, in years
    unitlength = 5000
    #unitlength = 2000
    #unitlength = 1000
    
    #J = 8 #silurian
    J = 18 #entire dataset
    #J = 24 #ordovician
    #J = 10 #LOME (450-448)

    nu = 4.24E-7*unitlength #entire dataset
    #nu = 3.28E-7*unitlength #ordovician
    #nu = 6.95E-7*unitlength #silurian
    #nu = 6.9E-7*unitlength #LOME

    s = 0.000495 #entire dataset
    #s = 0.0047 #silurian
    #s = 0.00021 #ordovician
    #s = 0.0001 #LOME

    timesteps = int(J*1E6/unitlength*20)
    trans = transmat(J=J,nu=nu)
    
    #empirical graptolite data
    datname = '' #replace with name of dataset with columns FA_age, LA_age representing first and last appearance in Ma
    df = pd.read_csv(datname)
    FADs = df.FA_age #you can also use df['column_name']
    LADs = df.LA_age
    ord = np.logical_and(list(df.FA_age>447),list(df.FA_age<481))
    ord_no_lome = np.logical_and(list(df.FA_age>449),list(df.FA_age<481))
    lome = np.logical_and(list(df.FA_age>448),list(df.FA_age<450))
    sil = np.logical_and(list(df.FA_age>419),list(df.FA_age<447))
    subset_old = 445
    subset_young = 440
    subset = np.logical_and(list(df.FA_age>subset_young),list(df.FA_age<subset_old))
    dur = (FADs-LADs)*1E6 #durations in years
    dur = dur/unitlength #duration in lifespans
    dur_ord = dur[ord]
    dur_sil = dur[sil]
    dur_ord_no_ME = dur[ord_no_lome]
    dur_ME = dur[lome]
    dur_subset = dur[subset]
    dur = dur_sil #uncomment and change this to only study subset
    dur = [math.ceil(i*J) for i in dur if i > 0] #this is durations in number of NT timesteps, rounded and removed of 0 durations
    #dur = [math.ceil(i) for i in dur if i > 0] #same, but in lifespans for comparing between different J
    se = surv(dur) #survivorship empirical
    se = np.divide(se,max(se))
    print(len(dur))

    #"""#theoretical
    start = time.time()
    haz,os = hazardandomegastar(tm=trans,Ni=1,t=timesteps,s=s) 
    print(time.time()-start,'seconds')
    #"""

    #"""
    #calculating log likelihoods for the whole spread
    start = time.time()
    delt = list(range(1,max(dur)+1)) #uncomment to plot entire duration
    timeinmy = 1
    #delt = list(range(1,int(J*1E6/unitlength*timeinmy))) #uncomment to plot part of duration
    ntprob = [0]+probdelt(delt=delt,d=os,h=haz,endpoint=timesteps,normalize=True)
    cdf = np.cumsum(ntprob)
    S = 1-cdf
    quantiles = [stats.binom(n=len(dur),p=float(i)).ppf(q=[0.025,0.975]) for i in S[1:]] #omit first entry because it's just 1
    lowerbound = [i[0]/len(dur) for i in quantiles]
    upperbound = [i[1]/len(dur) for i in quantiles]
    print(time.time()-start,'seconds')
    #print(ntprob)
    #print('first few probabilities',ntprob[:10])
    #print('ntprob total',sum(ntprob))
    print('log likelihood',sum([math.log(ntprob[i]) for i in dur]))   
    #"""
    
    """
    #calculating log likelihoods of only the observed durations
    delt = list(np.unique(dur))
    counts = [dur.count(i) for i in delt]
    start = time.time() 
    ntprob = probdelt(delt=delt,d=os,h=haz,endpoint=timesteps,normalize=True)
    print(time.time()-start,'seconds')
    print('log likelihood',sum([math.log(ntprob[i])*counts[i] for i in range(len(delt))]))
    #"""

    #"""
    #kolmogorov-smirnov test
    sorteddurations = np.unique(dur)
    ecdf = np.cumsum([sum([j==i for j in dur])/len(dur) for i in sorteddurations])
    tcdf = [cdf[i] for i in sorteddurations] #theoretical
    #plt.scatter(ecdf,tcdf)
    ks = max([abs(i) for i in np.subtract(tcdf,ecdf)])
    print('N =',len(dur))
    print('Kolmogorov-Smirnov stat',ks)
    print('Critical value, alpha = 0.05:',1.36/math.sqrt(len(dur)))
    print('Critical value, alpha = 0.01:',1.63/math.sqrt(len(dur)))
    #"""
    
    #"""
    #plotting
    fig = plt.figure()
    ax = plt.axes()
    plt.semilogy()
    plt.xlabel("Lifespans") 
    ax.plot([i/J for i in range(max(dur))],se)
    #"""

    #"""
    #plotting theoretical
    ax.plot([0]+[i/J for i in delt],S)
    #95% prediction interval
    ax.plot([i/J for i in delt],lowerbound)
    ax.plot([i/J for i in delt],upperbound)
    #"""

    #"""
    #weibull
    shape,scale = 0.5941937,4363.4815149 #entire dataset
    #shape,scale = 0.6082754,8170.766146 #ordovician
    #shape,scale = 0.6184011,1346.356428 #silurian
    #shape,scale = 1.401756,3181.640063 #LOME
    surv_weibull = [math.exp(-(i/scale)**shape) for i in range(max(dur))]
    quantiles_weibull = [stats.binom(n=len(dur),p=float(i)).ppf(q=[0.025,0.975]) for i in surv_weibull[1:]]
    lowerbound_weibull = [i[0]/len(dur) for i in quantiles_weibull]
    upperbound_weibull = [i[1]/len(dur) for i in quantiles_weibull]

    #exponential function
    scale_exp = 5939.265 #entire dataset
    #scale_exp = 10840.16 #ordovician
    #scale_exp = 1830.487 #silurian
    #scale_exp = 3479.546 #LOME
    surv_exp = [math.exp(-i/scale_exp) for i in range(max(dur))]
    quantiles_exp = [stats.binom(n=len(dur),p=float(i)).ppf(q=[0.025,0.975]) for i in surv_exp[1:]]
    lowerbound_exp = [i[0]/len(dur) for i in quantiles_exp]
    upperbound_exp = [i[1]/len(dur) for i in quantiles_exp]

    #plotting
    ax.plot([i/J for i in range(max(dur))],surv_weibull)
    ax.plot([i/J for i in range(1,max(dur))],lowerbound_weibull)
    ax.plot([i/J for i in range(1,max(dur))],upperbound_weibull)
    ax.plot([i/J for i in range(max(dur))],surv_exp)
    ax.plot([i/J for i in range(1,max(dur))],lowerbound_exp)
    ax.plot([i/J for i in range(1,max(dur))],upperbound_exp)
    #"""

    plt.show()