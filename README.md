# nt-surv
Python code for modeling age-dependent extinction with ecological neutral theory

markov.py file contains all the functions for calculating probabilities of extinction and sampling as a function of species age in Hubbell's neutral theory. This corresponds to equations 1-7 in the methods.

samplednt.py implements survivorship in neutral theory with incomplete sampling. This corresponds to equations 8-15 in the methods.

main.py is used for all the main analyses except fitting NT models, which is done in cumulated.py.

cumulated.py calculated log likelihoods for durations in units of lifespans, not timesteps. One NT individual has a lifespan equal to J timesteps, where J is community size. Therefore likelihoods cannot be compared for different J when data are in units of timesteps. This script rounds all data up to the nearest integer lifespan, and calculates probabilities accordingly. For example, for J=18, durations of between 1 and 18 timesteps get rounded up to 1 lifespan. So the corresponding probability for 1 lifespan cumulates (adds together) probabilities for between 1 and 18 timesteps.
