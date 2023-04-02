#markov.py file contains all the functions for calculating probabilities of extinction and sampling as a function of species age in Hubbell's neutral theory
#corresponds to equations 1-7 in the methods
import numpy as np, matplotlib.pyplot as plt, random, math, time
from scipy import linalg, sparse
from scipy.sparse.linalg import eigs
from matplotlib.patches import Rectangle

#helper function for transmat that returns indices of the kth off-diagonal for an a by a matrix 
def kth_diag_indices(a, k):
    rows, cols = np.diag_indices(a)
    if k < 0:
        return rows[-k:], cols[:k]
    elif k > 0:
        return rows[:-k], cols[k:]
    else:
        return rows, cols

#this function makes a markov matrix of size J+1 x J+1
#entries correspond to transition probabilty for a single local community with speciation and no migration
def transmat(J,nu):
    diag = np.diag_indices(J+1)
    diagr = kth_diag_indices(J+1,1)
    diagl = kth_diag_indices(J+1,-1)
    mat = np.zeros([J+1,J+1]) 
    #diagonal
    mat[diag] = \
    [(1-Ni/J)*((J-Ni-1)/(J-1) + Ni*nu/(J-1)) + Ni/J*(Ni-1)/(J-1)*(1-nu) for Ni in range(J+1)]
    #first right off-diag
    mat[diagr] = \
    [(J-Ni)/J*Ni/(J-1)*(1-nu) for Ni in range(J)]
    #first left off-diag
    mat[diagl] = \
    [Ni/J*((J-Ni)/(J-1)+(Ni-1)*nu/(J-1)) for Ni in range(1,J+1)]
    return(mat)

#get the probability of surviving to a certain age
def survivorship(tm,Ni,length):
    out = [1]
    tmm = np.identity(len(tm))
    for i in range(length):
        #tmm = np.linalg.matrix_power(tm,i)
        tmm = np.dot(tmm,tm)
        out.append(1-tmm[Ni][0])
    return(out)

#inducedmeasure() gives the induced measure, or probability distribution of abundance states
#terminology follows Kemeny & Snell 1960, Finite Markov Chains
#takes a transition matrix, a starting abundance, and the number of timesteps
#so this is the probability of a species that starts at abundance Ni being at each abundance after t timesteps
#NOTE: this and all its dependent functions have weird behavior for low values of J (below 10 or so) related to matrix_power
def inducedmeasure(tm,Ni,t):
    return(np.linalg.matrix_power(tm,t)[Ni])

#gives the induced measure (π*_Ni,t) after t timesteps conditional on the species still being extant, so P(N=0)=0
#this is so inefficient
def inducedmeasureextant(tm,Ni,t):
    im = inducedmeasure(tm,Ni,t)
    return(np.insert(im[1:]/(1-im[0]),0,0))

#scd() gives the stationary conditional distribution given a transition matrix
#from Darroch & Seneta 1965
def scd(tm):
    Q = tm[1:,1:]
    eval,evec = linalg.eig(Q,left=True,right=False)
    idx = eval.argsort()[::-1]   
    #eval = eval[idx]
    evec = evec[:,idx]
    return(evec[:,0]/sum(evec[:,0]))

#SAMPLED NEUTRAL THEORY
#omegastar gives the probability of observing a species extant at duration d at every interval of its existence t
#s is sampling probability per individual per timestep
def omegastar(tm,Ni,d,s):
    out = []
    #for every entry ω*_[Ni,t] in ω*_Ni, for t from 0 to d
    for t in range(d+1):
        pistar = inducedmeasureextant(tm=tm,Ni=Ni,t=t) #this is π*_[Ni,t]
        #prob. being at some abundance times prob of being observed given that abundance, summed over all possible abundances
        out.append(sum([pistar[x]*(1-(1-s)**x) for x in range(len(pistar))]))
    return(out)

#omega is like omegastar, but NOT conditional on survival
def omega(tm,Ni,d,s):
    out = []
    #for every entry ω*_Ni,t in ω*_Ni, t from 0 to d
    for t in range(d+1):
        pistar = inducedmeasure(tm=tm,Ni=Ni,t=t) #this is π*_Ni,t
        #prob. being at some abundance times prob of being observed given that abundance, summed over all possible abundances
        out.append(sum([pistar[x]*(1-(1-s)**x) for x in range(len(pistar))]))
    return(out)

#hazard function gives the probability of going extinct in the next timestep for each timestep P^t_1=>0
def hazardfunc(tm,Ni,t,every=1):
    return([np.dot(inducedmeasureextant(tm=tm,Ni=Ni,t=i),tm)[0] for i in range(0,t,every)])

#hazardandomegastar is a much faster way of jointly calculating hazard rate and omegastar
#returns a 2-tuple, first is hazard rate, second is omegastar
def hazardandomegastar(tm,Ni,t,s,every=1):
    haz=[]
    os=[]
    tmexp = np.identity(len(tm))
    for i in range(t):
        im = tmexp[Ni] #induced measure at time i
        cim = np.insert(im[1:]/(1-im[0]),0,0) #conditional induced measure
        haz.append(np.dot(cim,tm)[0])
        os.append(sum([cim[x]*(1-(1-s)**x) for x in range(len(cim))]))
        tmexp = np.dot(tmexp,tm) #P^i
    return((haz,os))

