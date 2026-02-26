Development Repository
================

## This is a development repository

It forks the [USE package](https://github.com/danddr/USE) with the goal to extend it to enable MCMC sampling in higher
dimensional spaces.

The repo can be downloaded by using the following command in a terminal or downloading the repo as a zip file:
```
git clone https://github.com/Mdmdma/USE.MCMC.git
```
The vignette [Insights on Markov chain based pseudo-absence sampling](https://mdmdma.github.io/USE.MCMC/articles/insights-on-MCMC-pseudo-absence-sampling-vignette.html) gives insights into the developed method, the vignette [Insights on nearest neighbor search](https://mdmdma.github.io/USE.MCMC/articles/Insights-on-nearest-neighbor-search.html) adds an explainer to the concept of using distance in the principal component space as a way to map fictional simulated points to real points in the environment. Both vignettes can be found on the github pages site.

#TLDR
We can use markov chains to sample pseudo absences in moderate dimensions. In the process we overcome limitations of the original paper, especially those related to a sampling grid, only two dimensions in the PCA space, and the hard cutoff threshold to exclude the presence species. 
