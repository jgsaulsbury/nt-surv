# nt-surv
Python code for modeling age-dependent extinction with ecological neutral theory

**markov.py** contains all the functions for calculating probabilities of extinction and sampling as a function of species age in Hubbell's neutral theory. This corresponds to equations 1-7 in the methods.

**samplednt.py** implements survivorship in neutral theory with incomplete sampling. This corresponds to equations 8-15 in the methods.

**main.py** is used for all the main analyses except fitting NT models, which is done in cumulated.py.

**cumulated.py** calculates log likelihoods for durations in units of lifespans, not timesteps. One NT individual has a lifespan equal to J timesteps, where J is community size. Therefore likelihoods cannot be compared for different J when data are in units of timesteps. This script rounds all data up to the nearest integer lifespan, and calculates probabilities accordingly. For example, for J=18, durations of between 1 and 18 timesteps get rounded up to 1 lifespan. So the corresponding probability for 1 lifespan cumulates (adds together) probabilities for between 1 and 18 timesteps.

**ntsim_local.py** and **ntsim_meta.py** provide functions for simulating neutral communities at the local community and metacommunity scale, respectively. They were not used in the paper, but are useful for exploring the predictions of neutral theory.

**supp_exercises.py** contains two analyses described in the supplementary materials: experiments with temporal resolution and with migration.

**supp_exercises.R** contains more experiments with temporal resolution.

**cramptonetal2016data.csv** is the dataset, used by main.py and cumulated.py.
