#function for simulating neutral community dynamics at different scales plus helpers and some examples
import random, matplotlib.pyplot as plt, numpy as np
from collections import Counter

"""
Simulating metacommunities
"""

#local communities in a dynamic metacommunity
#m1 and m2 are the dimensions of the metacommunity, which is a grid of local communities of size J
#m is migration rate, nu is speciation rate, t is number of timesteps to simulate for
#spcounter is used so that new species don't share a number with existing ones
#uses a few helper functions (see below)
def meta(m1,m2,J,m,nu,t,save_at=None):
    if save_at is None:
        save_at = [t] #returns only the last timestep by default
    meta_new = np.arange(m1*m2*J).reshape(m1,m2,J) #instantiates metacommunity with all singletons
    out = [meta_new] if save_at.count(0) > 0 else []
    spcounter = m1*m2*J #keeps track of how many species have been instantiated
    for iter in range(t): #iterates t times
        meta_old = meta_new
        meta_disturbed = disturb_meta(meta_old,m1,m2,J,1) #disturbs metacommunity
        meta_new = replace(meta_disturbed,m1,m2,J,m) #replaces with locals or migrants
        #convert additions to new species with probability nu
        for i in range(len(meta_new)):
            for j in range(len(meta_new[0])):
                if random.random() < nu: #if there is speciation
                    meta_new[i][j][J-1] = spcounter #set the new sp in this local community to be a new one
                    spcounter+=1
        if save_at.count(iter+1) > 0:
            out.append(meta_new)
    return(out)

#Helper functions for meta()

#flatten is a helper function for Moore()
def flatten(items, seqtypes=(list, tuple)):
    for i, x in enumerate(items):
        while i < len(items) and isinstance(items[i], seqtypes):
            items[i:i+1] = items[i]
    return items

#Moore() returns the Moore neighborhood of a local community
#given a >=2D matrix and a set of indices i and j, Moore() returns the concatenated contents of the indicated cell's Moore neighborhood 
#first gets the boundaries of the neighborhood
#slices those indices out, removes entry at i,j, then returns
def Moore(mat,i,j):
    top = max(0,i-1)
    bottom = min(len(mat),i+1)
    left = max(0,j-1)
    right = min(len(mat[0]),j+1)
    out = mat[top:bottom+1,left:right+1].tolist()
    out[i-top].pop(j-left)
    return(flatten(out))

#disturb() removes D individuals from a local community
#takes a list of species occurrences loc and disturbance size D
#so if loc is [1,1,1,2,1,4,4], that indicates 4 individuals of sp 1, 1 of sp 2, and 2 of sp 4
def disturb(loc,D):
    return(random.sample(list(loc),len(loc)-D))

#disturb_meta() takes a metacommunity, its dimensions and a disturbance size D and disturbs each local community
def disturb_meta(meta,m1,m2,J,D):
    return(np.array([disturb(meta[i][j],D) for i in range(m1) for j in range(m2)]).reshape(m1,m2,J-D))

#replace takes a disturbed metacommunity and its undisturbed dimensions
#returns a "filled" metacommunity with only local replacement
#each cell is replaced by migrants from the Moore neighborhood with probability m
#for simplicity of code D is set to 1
#if-else statement returns a migrant with probability m, else does what replace_local does
def replace(meta,m1,m2,J,m):
    return(np.array([[*meta[i][j],*random.choices(Moore(meta,i,j),k=1)] if random.random() < m else \
        [*meta[i][j],*random.choices(meta[i][j],k=1)] \
        for i in range(m1) for j in range(m2)]).reshape(m1,m2,J))


#Functions for visualizing the output of meta()

#map_species visualizes a species's abundances across the map
def map_species(meta,sp):
    out = np.array([sum(meta[i][j] == sp) for i in range(len(meta)) for j in range(len(meta[0]))]).reshape(len(meta),len(meta[0]))
    return(out)

#count_species
def count_species(meta):
    return(len(np.unique(meta)))

#trajectory() takes a list of metapopulations and a species ID and outputs the number of individuals of that species at each timestep
def trajectory(metalist,sp):
    return([np.count_nonzero(meta==sp) for meta in metalist])

#rankabundance returns an ordered (large-small) list of the abundances of all the sp in a metacommunity
def rankabundance(meta):
    out = list(Counter(meta.flatten()).values())
    out.sort(reverse=True)
    return(out)

"""
Main loop with some examples
"""
if __name__ == '__main__':
    random.seed(1)

    x = meta(m1=8,m2=8,J=10,m=0.05,nu=0.01,t=6000,save_at=range(0,6001,1)) #low migration
    #print(x[6000])
    print(map_species(x[6000],3271)) #where in the 8x8 metacommunity is species 3271
    print(count_species(x[6000]))
    print(rankabundance(x[6000]))
    plt.plot(trajectory(x,3271)) #plot total abundance of species 3271 through time
    plt.plot(trajectory(x,2001)) #with an extinct one for comparison
    plt.show()