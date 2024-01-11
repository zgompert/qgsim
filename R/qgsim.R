qgsim <-
function(npops=2,mig=0,Ne=500,theta0=0,z0=0,h2=0.5,Gcor=0.2,omega11=1,omega22=1,omegaCor=0,model="Brownian",
    ngens=100,tsd=0.02,tmn=0.01){

	## always consider two traits
	ntraits<-2    
    
	## construct the Gmatrix from the trait heritabilities, assumed to be the same for both traits,
	## and genetic correlations
	## this assumes the phenotypic variance is 1
	G<-matrix(c(h2,Gcor,Gcor,h2),nrow=2,byrow=TRUE)
	## construct omega matrix, curvature of the adaptive landscape
	omega<-matrix(c(omega11,omegaCor,omegaCor,omega22),nrow=2,byrow=TRUE)
    
	## users can give one value, one value per trait, or a full matrix of thetas
	## this creates a matrix from whatever they give
	theta0<-matrix(theta0,nrow=npops,ncol=ntraits,byrow=TRUE)
	## same for initial trait values, z0
	z0<-matrix(z0,nrow=npops,ncol=ntraits,byrow=TRUE)
    
	## create objects to store trait means and optimal trait values over time
	## one element in a list per population, each element is a matrix
	## rows = generations, columns = 2 traits
	z_mat<-vector("list",npops)
	theta_mat<-vector("list",npops)
    
	## initialize the z and theta matrixes with initial values from z0 and theta0
	for(j in 1:npops){
		z_mat[[j]]<-matrix(z0[j,],nrow=ngens,ncol=ntraits)
		theta_mat[[j]]<-matrix(theta0[j,],nrow=ngens,ncol=ntraits)
	}

	## evolution, loop over generations and popualtions
	for(i in 2:ngens){
		for(j in 1:npops){
            		## for conveince, grab current values of z and theta
			z<-z_mat[[j]][i-1,]
			theta<-theta_mat[[j]][i-1,]
            		## update theta, optimal trait values, based on the dynamic
			## adaptive landscape model chosen by the user
			if(model=="Brownian"){ ## Brownian random walk
		                theta<-theta+rnorm(ntraits,0,tsd)
           		} else if (model=="Uncorrelated"){ ## sample normal each time, independent
				theta<-rnorm(traits,0,tsd)
			} else if (model=="Trend"){ ## changes by trend plus noise
				theta<-rnorm(traits,tmn,tsd)
			}
			## evolution happens, selection and drift	
			## z' = z + (theta - z)/omega * G + N(0,G/Ne)
			## matrix algebra fun
			z<-z + (theta-z) %*% solve(omega) %*% G + rmnorm(1,varcov=G/Ne)
			## store values
			z_mat[[j]][i,]<-z
			theta_mat[[j]][i,]<-theta
		}
		## add gene flow, which happens after selection in this model
		z1_vec<-rep(NA,npops)
		z2_vec<-rep(NA,npops)
		for(j in 1:npops){
			z1_vec[j]<-z_mat[[j]][i,1]
			z2_vec[j]<-z_mat[[j]][i,2]
		}
		z1_bar<-mean(z1_vec);z2_bar<-mean(z2_vec)
		## z' = z (1-m) + z_bar * m
		for(j in 1:npops){
			z_mat[[j]][i,1]<-z_mat[[j]][i,1] * (1-mig) + z1_bar * mig
			z_mat[[j]][i,2]<-z_mat[[j]][i,2] * (1-mig) + z2_bar * mig
		}

    }
    o<-list(z=z_mat,theta=theta_mat)
    return(o)
}
