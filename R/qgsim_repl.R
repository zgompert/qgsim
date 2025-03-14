#' Simulate evolution (with replicates and summaries)
#'
#' This function is a wrapper for the gqsim function which lets you run replicates and get summaries.
#' The summaries currently available include: mean rate of evolution ('mean_evol_rate', average change in each trait across populations and generations), evolutionary lag ('mean_evol_lag', average distance in bivaraiate trait space between the optimal phenotypes and z), correlations in evolution between traits ('mean_traitwise_corr', averaged over time and populations), correlations in evolution across populations ('mean_popwise_corr', averaged over time and traits), and the variance in mean trait values across populations in the final generation ('end_phenotypic_variance', averagec over traits).
#' @param nreps number of replicates to run
#' @param summaries a vector of strings containing any summaries to be included in the results
#' @param npops number of populations to simulate.
#' @param mig the migration rate, this is the total proportion of individuals in each population made up of migrants from the other populations.
#' @param theta0 this specifies the initial optimal value for the traits, you can supply a single value (applied to both traits), one value per trait, or a matrix with one value per population (row) and trait (column).
#' @param z0 this specifies the initial trait values, you can supply a single value (applied to both traits), one value per trait, or a matrix with one value per population (row) and trait (column).
#' @param h2 the trait heritabilities, assumed to be the same for both traits (must be between 0 and 1).
#' @param Gcor the genetic correlation between the pair of traits (must be between -1 and 1).
#' @param omega11 denotes the intensity of selection (curvature of the adaptive landscape) with respect to trait 1, larger values result in weaker selection [default = 1].
#' @param omega22 denotes the intensity of selection (curvature of the adaptive landscape) with respect to trait 2, larger values result in weaker selection [default = 1].
#' @param omegaCor denotes the strength of correlational selection, larger values denote stronger selection for combinations of trait 1 and trait 2 (must be between -1 and 1, 0 denotes independent selection on each trait) [default = 0].
#' @param model denotes the model for adaptive peak movement, must be one of the following: "Brownian", "Uncorrelated", or "Trend"; detailed descriptions are provided below.
#' @param ngens number of generations to simulate [default = 100].
#' @param tsd standard deviation for peak movement, larger values denote larger (random) jumps in adaptive peak locations (must be a positive number).
#' @param tmn average direction shift in the location of the adaptive landscape, only relevant for the "Trend" model (can be negative or positive, but must be a single value for both traits).
#'
#' @details
#' With the Brownian motion peak shift model, a random normal deviation with standard deviation tsd is added to the peak location each generation. With the Uncorrelated peak shift model, the peak location is follows a normal distribution with independent peak locations each generation. The peak location for each trait is specifically drawn from a normal distribution with mean 0 and standard deviation tst. Finally for Trend, the peak location shifts by a random normal deviation each generation with a mean tmn and standard deviation tsd. In all cases, peak shifts occur independently in each population.
#'
#' @return A list with two objects, the mean trait values (z) and adaptive peak locations (theta). Each of these is itself a list with one matrix per population. The population matrixes have two columns, one per trait, and one row per generation.
qgsim_repl <- function(nreps=1, summaries=c(), npops=2,mig=0,Ne=500,theta0=0,z0=0,h2=0.5,Gcor=0.2,
    omega11=1,omega22=1,omegaCor=0,model="Brownian",ngens=100,tsd=0.02,tmn=0.01) {

  # check that requested summaries are valid
  valid_summaries <- c('mean_evol_rate', 'mean_evol_lag', 'mean_traitwise_corr', 'mean_popwise_corr', 'end_phenotypic_variance')
  invalid_summaries <- setdiff(summaries, valid_summaries)
  if (length(invalid_summaries) > 0) {
    stop(paste("Error: Invalid summary/summaries: ", paste(invalid_summaries, collapse=", "), ". Must be in: ", paste(valid_summaries, collapse = ", "), ".", sep=""))
  }

  res.z<-list()
  res.theta<-list()
  
  ncols = length(summaries)+1
  sum_mat <- matrix(nrow=nreps, ncol=ncols)
  colnames(sum_mat) <- c('replicate', summaries)
  sum_mat[,'replicate']<-1:nreps
  
  for (i in 1:nreps) {
    res<-qgsim_func(npops=npops,mig=mig,Ne=Ne,theta0=theta0,z0=z0,h2=h2,Gcor=Gcor,omega11=omega11,omega22=omega22,omegaCor=omegaCor,model=model,ngens=ngens,tsd=tsd, tmn=tmn)
    res.z[[i]]<-res$z
    res.theta[[i]]<-res$theta
  }

  # add in summaries

  if ("mean_evol_rate" %in% summaries) {
    # average rate of change in each trait across generations and populations

    # sum the differences in each trait
    for (rep_i in 1:nreps) {
      delta_sum <- 0
      z <- res.z[[rep_i]]
      for (pop_i in 1:npops) {
        pop <- z[[pop_i]]
        for (gen_i in 2:ngens) {
          # add trait 1 delta
          delta_sum <- delta_sum + pop[gen_i, 1] - pop[gen_i-1, 1]
          # add trait 2 delta
          delta_sum <- delta_sum + pop[gen_i, 2] - pop[gen_i-1, 2]
        }
      }
      # divide by (# traits) * (# generation deltas) * (# populations)
      mean_evol_rate <- delta_sum / (2 * (ngens-1) * npops)
      sum_mat[rep_i,'mean_evol_rate'] <- mean_evol_rate
    }
  }
  
  if ("mean_evol_lag" %in% summaries) {
    for (rep_i in 1:nreps) {
      delta_sum <- 0
      z <- res.z[[rep_i]]
      theta <- res.theta[[rep_i]]
      for (pop_i in 1:npops) {
        pop.z <- z[[pop_i]]
        pop.theta <- theta[[pop_i]]

        t1_dist <- sum(abs(pop.z[,1] - pop.theta[,1]))
        t2_dist <- sum(abs(pop.z[,2] - pop.theta[,2]))
        delta_sum <- delta_sum + t1_dist + t2_dist
      }
      # divide by (#traits) * (# generation deltas) * (# populations)
      mean_evol_lag <- delta_sum / (2 * ngens * npops)
      sum_mat[rep_i, 'mean_evol_lag'] <- mean_evol_lag
    }
  }

  if ("mean_traitwise_corr" %in% summaries) {
    for (rep_i in 1:nreps) {
      corr_sum <- 0
      z <- res.z[[rep_i]]
      for (pop_i in 1:npops) {
        pop.z <- z[[pop_i]]
        # adding correlation between trait values in a population
        corr_sum <- corr_sum + cor(pop.z[,1], pop.z[,2])
      }
      # divide by (# populations)
      mean_traitwise_corr <- corr_sum / npops
      sum_mat[rep_i, 'mean_traitwise_corr'] <- mean_traitwise_corr
    }
  }

  if ("mean_popwise_corr" %in% summaries) {
    for (rep_i in 1:nreps) {
      corr_sum <- 0
      ncors<-0
      z <- res.z[[rep_i]]
      for (pop1_i in 1:(npops-1)) {
        for (pop2_i in (pop1_i+1):npops) {
          pop1.z <- z[[pop1_i]]
          pop2.z <- z[[pop2_i]]
          corr_sum <- corr_sum + cor(pop1.z[,1], pop2.z[,1])
          corr_sum <- corr_sum + cor(pop1.z[,2], pop2.z[,2])
          ncors <- ncors + 2
        }
      }
      # divide by number of correlations
      mean_popwise_corr <- corr_sum / ncors
      sum_mat[rep_i, 'mean_popwise_corr'] <- mean_popwise_corr
    }
  }

  if ("end_phenotypic_variance" %in% summaries) {
    for (rep_i in 1:nreps) {
      z <- res.z[[rep_i]]
      z_end<-matrix(NA, nrow=npops, ncol=2) ## two traits
      for (pop in 1:npops) {
          z_end[pop,1]<-z[[pop]][ngens,1]
          z_end[pop,2]<-z[[pop]][ngens,2]      
      }
      # calcualte variance and take mean for the two traits
      end_pheno_var <-mean(apply(z_end,2,var))
      sum_mat[rep_i, 'end_phenotypic_variance'] <- end_pheno_var
    }
  }

  return(sum_mat)
}

