#Simulation of species presence and presence-absence data for ecological niche modeling (ENM), also known as species distribution modeling (SDM).

 Ecological niche models (ENM) or species distribution models (SDM) attempt to predict the spatiotemporal distribution of species on the basis of environmental data and species observation records In addition to simply delineating the distribution of a species, ENMs often output an index representing either local relative suitability of habitat (using presence-only observations) or local probabilities of species presence (using presence-absence observations). The distinction between presence-only and presence-absence species observation records is important because it impacts options for analysis and as described above, downstream interpretation of results.

 The aim of this tutorial is to simulate environmentally determined species presence and presence-absence data in R. Specifically, this tutorial covers how to 1) simulate large numbers of realistic environmental surfaces, 2) simulate a true species presence probability surface from those environmental surfaces, and 3) simulate species presence or presence-absence records from the true species presence probability surface. Typically, the goal of ENMs/SDMs is to recover both the true species presence probability surface and understand the relationships between that surface and environmental data. Accordingly, as the truth is known, the outputs of this simulation pipeline can be used to test the efficacy of ENM methods.