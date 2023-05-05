#cumulated.py calculated log likelihoods for durations in units of lifespans, not timesteps
#one NT individual has a lifespan equal to J timesteps, where J is community size
#therefore likelihoods cannot be compared for different J when data are in units of timesteps
#Rounds all data up to the nearest integer lifespan, and calculates probabilities accordingly
#For example, for J=18, durations of between 1 and 18 timesteps get rounded up to 1 lifespan
#So the corresponding probability for 1 lifespan cumulates (adds together) probabilities for between 1 and 18 timesteps

import time, numpy as np, math, pandas as pd
from markov import transmat, hazardandomegastar, probcumulate
from samplednt import probdelt
from main import surv

#main loop
if __name__ == '__main__':
    #params
    unitlength = 5000 #a graptolite unit has this expected lifespan
    J = [17,18,19] #all the community sizes to consider
    s = [0.000495,0.00049,0.0005] #all the sampling rates to consider
    #ordovician
    #nu = 4.24E-7*unitlength
    #silurian
    #nu = 6.95E-7*unitlength
    #ordovician minus LOME (481-449)
    #nu = 3.05E-7*unitlength
    #entire dataset
    nu = 4.24E-7*unitlength
    #other
    #nu = 6.44E-7*unitlength

    #empirical graptolite data
    datname = 'cramptonetal2016data.csv'
    df = pd.read_csv(datname)
    FADs = df.FA_age #you can also use df['column_name']
    LADs = df.LA_age
    ord = np.logical_and(list(df.FA_age>447),list(df.FA_age<481))
    lome = np.logical_and(list(df.FA_age>448),list(df.FA_age<450))
    sil = np.logical_and(list(df.FA_age>419),list(df.FA_age<447))
    dur = (FADs-LADs)*1E6 #durations in years
    #dur = dur[ord] #uncomment for only ordovician
    #dur = dur[sil] #uncomment for only silurian
    #dur = dur[lome] #uncomment for only LOME
    dur = dur/unitlength #duration in lifespans
    dur = [math.ceil(i) for i in dur if i > 0] #duration in lifespans (rounded) and removed of 0 durations
    se = surv(dur) #survivorship empirical
    se = np.divide(se,max(se))

    #log likelihoods
    loglik = pd.DataFrame(index=['J'+str(round(i,5)) for i in J],columns=['s'+str(round(i,5)) for i in s])  #store likelihoods here
    print(loglik)
    print('---')
    estimated = np.nan
    for x in range(len(J)):
        trans = transmat(J=J[x],nu=nu)
        for y in range(len(s)):
            timesteps = int(J[x]*1E6/unitlength*20) #calculate to 20 million years
            haz,os = hazardandomegastar(tm=trans,Ni=1,t=timesteps,s=s[y]) 
            print('J:',J[x],'s:',s[y])
            start = time.time() 
            delt = list(np.unique(dur))
            delt_timesteps = [item for sublist in [list(range(i*J[x]-J[x]+1,i*J[x]+1)) for i in delt] for item in sublist]
            ntprob = probdelt(delt=delt_timesteps,d=os,h=haz,endpoint=timesteps,normalize=True)
            if not np.isnan(float(ntprob[0])):
                ntprob_lifespans = probcumulate(ntprob,divisor=J[x])
                counts = [dur.count(i) for i in delt]
                lik = sum([math.log(ntprob_lifespans[i])*counts[i] for i in range(len(delt))])
                print('log likelihood',lik)
                end = time.time()
                print(end-start,'seconds')
                print('estimated time remaining:',(end-start)*(len(J)*len(s) - x*len(s) - y - 1),'seconds')
                loglik.iloc[x,y] = lik
                print('---')
    print(loglik)
    loglik.to_csv('loglik.csv')
