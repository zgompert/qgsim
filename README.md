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
Run and summarize multiple replicate simulations with `qgsim_repl`.

```
## load the qgsim package
library(qgsim)
## type help(qgsim_repl) to view the help page

## ten replicates each with h2 = 0.01, 0.5 or 0.99
## compute and report the average lag between optimal and mean trait values
qgsim_repl(nreps=10,c("mean_evol_lag"),npops=2, mig=0, Ne=500, h2=0.01, Gcor=0, omega11=1, omega22=1, omegaCor=0, model="Trend", ngens=100, tmn=.05,tsd=0.05)
qgsim_repl(nreps=10,c("mean_evol_lag"),npops=2, mig=0, Ne=500, h2=0.5, Gcor=0, omega11=1, omega22=1, omegaCor=0, model="Trend", ngens=100, tmn=.05,tsd=0.05)
qgsim_repl(nreps=10,c("mean_evol_lag"),npops=2, mig=0, Ne=500, h2=0.99, Gcor=0, omega11=1, omega22=1, omegaCor=0, model="Trend", ngens=100, tmn=.05,tsd=0.05)

## see the help page for other valid summaries
```
Run a single simulation with `qgsim_func` and then look at plots of the optimal value versus the mean trait value for each trait with `qgsim_plot_bivar`. This function can plot pairs of z or theta for one or two populations.

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

## plot of optimal versus mean for trait 1
qgsim_plot_bivar(simOut,1,c("z1","theta1"),type='l')

## plot of optimal versus mean for trait 2
qgsim_plot_bivar(simOut,1,c("z2","theta2"),type='l')
```
Use the helper function `var_covar_mat_vis` to visualize G and Omega matrixes. This can be useful to make example figures or to aid in your decision making about parameters to use.

```R
## load the qgsim package
library(qgsim)
## make the plot
var_covar_mat_vis(h2=.5,Gcor=.3,omega11=1,omega22=1,omegaCor=-0.3)
```
