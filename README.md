# qgsim

# Installation

You can install this package within R directly from GitHub.

```R
## install and load devtools
install.packages("devtools")
library(devtools)
## install qgwim
devtools::install_github("zgompert/qgsim")
## load the qgsim package
library(qgsim)
```
# Usage

See the manual [qgsim_1.0.pdf](qgsim_1.0.pdf) or help pages for detailed descriptions of each function.

Brief examples illustrating usage of the core functions are provided below.

Run a single simulation with `qgsim_func`.

```R
## load the qgsim package
library(qgsim)
## type help(qgsim_func) to view the help page

## simulation with:
## 2 populations, no migration,
## h2 = .5, genetic correlation = 0
## omega11 = omega22 = 1, omegaCor = 0
## Brownnian model, tsd = .05
simOut<-qgsim_func(npops=2, mig=0, Ne=500, h2=0.5, Gcor=0, omega11=1, omega22=1, omegaCor=0, model="Brownian", ngens=100, tsd=0.05)

## view optimal trait values (theta) for population 1
simOut$theta[[1]]

## view mean trait values (z) for population 1
simOut$z[[1]]

## view lag, mean - optimal, for trait 1
lag<-simOut$z[[1]]-simOut$theta[[1]]
lag

## compute mean lag over simulation for each trait
apply(lag,MARGIN=2,mean)

```
