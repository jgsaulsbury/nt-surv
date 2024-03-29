#assorted supplemental exercises
from ntsim_meta import *
from markov import *
from main import *
import matplotlib.pyplot as plt, pandas as pd,mpmath

#main loop
if __name__ == '__main__':
        
    #remove pound sign on next line to comment exercise out
    #"""==EXERCISE 1: coarsening temporal resolution==
    unitlength = 5000
    J = 18 #entire dataset
    #nu = 4.56E-7*unitlength #entire dataset
    nu = 4.24E-7*unitlength #not rounded up
    #s = 0.00037 #entire dataset
    s = 0.000495 #not rounded up
    #empirical graptolite data
    datname = 'cramptonetal2016data.csv'
    df = pd.read_csv(datname)
    FADs = df.FA_age #you can also use df['column_name']
    LADs = df.LA_age
    dur = (FADs-LADs)*1E6 #durations in years
    dur = dur/unitlength #durations in lifespans
    dur = [math.ceil(i*J) for i in dur if i > 0] #durations in number of NT timesteps, rounded and removed of 0 durations
    roundto = 1000000 #round up to this resolution
    resolution = round(roundto/unitlength*J) #resolution in timesteps = resolution in years/unitlength*J
    dur_coarse = [math.ceil(i/resolution)*resolution for i in dur] #coarsened resolution
    se = surv(dur_coarse) #survivorship empirical
    se = np.divide(se,max(se))
    #theoretical
    timesteps = int(J*1E6/unitlength*20)
    trans = transmat(J=J,nu=nu)
    start = time.time() #start timer
    haz,os = hazardandomegastar(tm=trans,Ni=1,t=timesteps,s=s) 
    #plotting
    #theoretical
    start=time.time()
    print("here goes")
    delt_plot = list(range(1,max(dur)+1))
    ntprob_plot = [0]+probdelt(delt=delt_plot,d=os,h=haz,endpoint=timesteps,normalize=True)    
    S = 1-np.cumsum(ntprob_plot)
    print(time.time()-start)
    #original resolution
    x = [0]+[i/J*unitlength/1E6 for i in delt_plot]
    plt.semilogy(x,S)
    plt.plot([0,15,15,0,0],[1,1,1E-4,1E-4,1]) #box for aligning with R output
    plt.show()
    #0.05My
    res = round(50000/unitlength*J)
    plt.semilogy([i for i in x[0::res] for _ in (0,1)][1:]+[x[-1]],[i for i in S[0::res] for _ in (0,1)])
    plt.plot([0,15,15,0,0],[1,1,1E-4,1E-4,1])
    plt.show()
    #0.1My
    res = round(100000/unitlength*J)
    plt.semilogy([i for i in x[0::res] for _ in (0,1)][1:]+[x[-1]],[i for i in S[0::res] for _ in (0,1)])
    plt.plot([0,15,15,0,0],[1,1,1E-4,1E-4,1])
    plt.show()
    #0.25My
    res = round(250000/unitlength*J)
    plt.semilogy([i for i in x[0::res] for _ in (0,1)][1:]+[x[-1]],[i for i in S[0::res] for _ in (0,1)])
    plt.plot([0,15,15,0,0],[1,1,1E-4,1E-4,1])
    plt.show()
    #0.5My
    res = round(500000/unitlength*J)
    plt.semilogy([i for i in x[0::res] for _ in (0,1)][1:]+[x[-1]],[i for i in S[0::res] for _ in (0,1)])
    plt.plot([0,15,15,0,0],[1,1,1E-4,1E-4,1])
    plt.show()
    #1My
    res = round(1000000/unitlength*J)
    plt.semilogy([i for i in x[0::res] for _ in (0,1)][1:]+[x[-1]],[i for i in S[0::res] for _ in (0,1)])
    plt.plot([0,15,15,0,0],[1,1,1E-4,1E-4,1])
    plt.show()
    #"""

    """==EXERCISE 2: migration==
    #simulates NT in an 8x8 metacommunity with J=8 and varying migration rate m
    random.seed(1)
    m1=m2=8
    J=8
    nu=0.005
    t = 50000
    #no migration
    x1 = meta(m1=m1,m2=m2,J=J,m=0,nu=nu,t=t,save_at=range(0,t+1,1)) #x1 stores the metacommunity at each timestep between 0 and t
    durs1 = [duration(metalist=x1,sp=i) for i in range(m1*m2*J,np.amax(x1[t]))]
    su1 = survivorshipemp(durs1)
    su1 = [i/max(su1) for i in su1]
    #low migration
    x2 = meta(m1=m1,m2=m2,J=J,m=0.02,nu=nu,t=t,save_at=range(0,t+1,1))
    durs2 = [duration(metalist=x2,sp=i) for i in range(m1*m2*J,np.amax(x2[t]))]
    su2 = survivorshipemp(durs2)
    su2 = [i/max(su2) for i in su2]
    #medium migration
    x3 = meta(m1=m1,m2=m2,J=J,m=0.1,nu=nu,t=t,save_at=range(0,t+1,1))
    durs3 = [duration(metalist=x3,sp=i) for i in range(m1*m2*J,np.amax(x3[t]))]
    su3 = survivorshipemp(durs3)
    su3 = [i/max(su3) for i in su3]
    #high migration
    x4 = meta(m1=m1,m2=m2,J=J,m=0.5,nu=nu,t=t,save_at=range(0,t+1,1))
    durs4 = [duration(metalist=x4,sp=i) for i in range(m1*m2*J,np.amax(x4[t]))]
    su4 = survivorshipemp(durs4)
    su4 = [i/max(su4) for i in su4]
    plt.semilogy([i/J for i in range(800*8+1)],survivorship(tm=transmat(J=8,nu=nu),Ni=1,length=800*8))
    plt.semilogy([i/(J*m1*m2) for i in range(800*512+1)],survivorship(tm=transmat(J=512,nu=0.005),Ni=1,length=800*512))
    plt.semilogy([i/J for i in range(len(su1))],su1)
    plt.semilogy([i/J for i in range(len(su2))],su2)
    plt.semilogy([i/J for i in range(len(su3))],su3)
    plt.semilogy([i/J for i in range(len(su4))],su4)
    plt.show()
    #"""
