#this file aggregates two different functions for simulating neutral local community dynamics
#these approaches store communities in different ways, so best not to assume they are interoperable
import random, copy, numpy as np
from collections import Counter

"""
Simulating local communities
"""

#ntlocal() simulates a local community with no metacommunity
#takes a local community (tuple of lists) and advances it
#first entry in the tuple gives the identities of the species in the community
#second entry gives the abundance of those species
#this way of storing the community is faster for some applications
#nu is the speciation probability, nsp is the total number of species thusfar
#disturb, replace, speciate
def ntlocal(comm,J,nu,timesteps,every=1): 
    nsp = max(comm[0])
    out = []
    while timesteps > 0: 
        killindex = random.randint(1,J) #kill this individual
        #print('killindex',killindex)
        i=0 #iterate this index along the community until you reach the killindex
        while killindex > sum(comm[1][0:i+1]):
            i += 1
        if comm[1][i] == 1:#if that species has 1 individual left
            del comm[0][i]
            del comm[1][i] #it goes extinct
        else: #otherwise just decrement
            comm[1][i] -= 1 
        repindex = random.randint(1,J-1) #this individual reproduces
        #print('repindex',repindex)
        j = 0 #iterate this index along the community until you reach the repindex
        while repindex > sum(comm[1][0:j+1]):
            j += 1
        if random.uniform(0,1)>nu: #if no speciation
            comm[1][j] += 1 #increment reproducing species
        else: #if speciation
            nsp += 1 
            comm[0].append(nsp)
            comm[1].append(1)
        timesteps -= 1
        if(timesteps % every == 0):
            out.append(copy.deepcopy(comm))
    return(out)

#local community embedded in static metacommunity
#takes [comm] (a list of integers representing the abundance of each species) and advances it [timesteps] cycles
#this is a model of the abundance dynamics of a local community embedded in a static metacommunity
#the 0th index of comm represents a null species
#out returns the state of the community at [every] timestep
#m is migration rate, J is local community size
#[metacomm] is a list of relative abundances of different species
#if comm is None, creates a local community with only the null species#
def ntlocal_in_meta(m,J,metacomm,timesteps,every=1,comm=None):
    if comm is None:
        comm = ([0]*(len(metacomm)-1))
        comm.insert(0,J)
    metacomm_cumsum = np.cumsum(metacomm) #for determining migrants
    out = []
    current_timestep = 0
    while current_timestep < timesteps: #in every timestep, 
        killindex = random.randint(1,J) #select a random individual to kill
        i=0 #increment species index i along the community until you reach the killindex
        while killindex > sum(comm[0:i+1]):
            i += 1
        comm[i] -= 1 #decrement that species
        if random.uniform(0,1)>m: #if no migration
            #this individual reproduces
            repindex = random.randint(1+comm[0],J-1) if comm[0<J-1] else 0 #uncomment to prevent null species from reproducing
            #repindex = random.randint(1,J-1)
            j = 0 #increment species index along the community until you reach the repindex
            while repindex > sum(comm[0:j+1]):
                j += 1
            comm[j] += 1 #increment reproducing species
        else: #otherwise if there is migration
            migrateindex = random.uniform(0,1)#select a migrant from the metacommunity
            k=0 #increment species index i along the metacommunity until you reach the killindex
            while migrateindex > metacomm_cumsum[k]:
                k += 1
            comm[k] += 1 #increment that species in the local community
        current_timestep += 1
        if(current_timestep % every == 0):
            out.append(copy.deepcopy(comm))
    return(out)

"""
Some helper functions for ntlocal_in_meta
"""

#metacomm generates a metacommunity using Fisher's log-series
#constant should be between 0 and 1
def meta(div,constant):
    out = [(-constant**i)/(i*np.log(1-constant)) for i in range(1,div+1)]
    out = [j/sum(out) for j in out]
    out.insert(0,0)
    return(out)

#commdiv() returns the diversity of a local community, minus the null species
def commdiv(comm):
    return(sum([i>0 for i in comm[1:]]))


"""
Main loop with some examples
"""
if __name__ == '__main__':
    random.seed(1)
    
    #parameters
    J=16
    nu=0.1
    t = 100

    #ntlocal
    c0 = ([0],[J])
    x = ntlocal(comm=c0,J=J,nu=nu,timesteps=t)
    #print(x)

    #ntlocal_in_meta
    m = 0.1
    metacomm = meta(div=10,constant=0.95)
    y = ntlocal_in_meta(m=m,J=J,metacomm=metacomm,timesteps=t,every=1)
    print(y)

