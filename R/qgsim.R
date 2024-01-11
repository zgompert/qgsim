qgsim <-
function(npops=2,mig=0,Ne=500,theta0=0,z0=0,h2=0.5,Gcor=0.2,omega11=1,omega22=1,omegaCor=0,model="Brownian",
    ngens=100,tsd=.02){

    ntraits<-2    
    
    G<-matrix(c(h2,Gcor,Gcor,h2),nrow=2,byrow=TRUE)
    omega<-matrix(c(omega11,omegaCor,omegaCor,omega22),nrow=2,byrow=TRUE)
    
    theta0<-matrix(theta0,nrow=npops,ncol=ntraits,byrow=TRUE)
    z0<-matrix(z0,nrow=npops,ncol=ntraits,byrow=TRUE)
    
    z_mat<-vector("list",npops)
    theta_mat<-vector("list",npops)
    
    for(j in 1:npops){
        z_mat[[j]]<-matrix(z0[j,],nrow=ngens,ncol=ntraits)
        theta_mat[[j]]<-matrix(theta0[j,],nrow=ngens,ncol=ntraits)
    }

    for(i in 2:ngens){
        for(j in 1:npops){
            z<-z_mat[[j]][i-1,]
            theta<-theta_mat[[j]][i-1,]
            if(model=="Brownian"){
                theta<-theta+rnorm(ntraits,0,tsd)
            }
            ## z' = z + (theta - z)/omega * G + N(0,G/Ne)
            z<-z + (theta-z) %*% solve(omega) %*% G + rmnorm(1,varcov=G/Ne)
            z_mat[[j]][i,]<-z
            theta_mat[[j]][i,]<-theta
        }
    }
    o<-list(z=z_mat,theta=theta_mat)
    return(o)
}
