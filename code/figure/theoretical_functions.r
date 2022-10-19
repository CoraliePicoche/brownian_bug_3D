library(pracma)


# Thomas and Poisson processes

#Ripley's functions
thomas_cdf=function(r,sigma,conc){
        a=((4/3)*pi*r^3)+1/(sigma*sqrt(pi)*conc)*(sqrt(pi)*sigma*erf(r/(2*sigma))-r*exp(-(r^2/(4*sigma^2))))
        return(a)
}


## Brownian Bug Model functions

#Ripley's K-functions
BBM_cdf=function(r,gamma,Delta,Conc,lambda,tau,t=NA){
        D=(Delta^2)/(2*tau)
        res=NA
        A=4/3*pi*r^3
        t=t*tau
        if(gamma>0){
                B=r^2/(6*D)
                C=(sqrt(gamma)/(6*sqrt(2)*(D^1.5)))*(r^3)*atan(r*sqrt(gamma/(2*D)))
                Dd=log((gamma*r^2)/(2*D)+1)/(6*gamma)
                E=(sqrt(gamma)*pi*r^3)/(12*sqrt(2)*D^1.5)
                res=A+(2*lambda/Conc)*(B+C+Dd-E)
        }else{
                B=r^2/2
                C=1/2*erf(r/(sqrt(8*D*t)))*(r^2-4*D*t)
                Dd=sqrt(2*D*t)*r/sqrt(pi)*exp(-r^2/(8*D*t))
                res=A+(lambda/(Conc*D))*(B-C-Dd)
        }
        return(res)
}

#Pcf
G_theoretical=function(gamma,r,lambda,Delta,C_0,tau,Tmax=NA){ #U=0 in the absence of advection
        D=(Delta^2)/(2*tau)
        Tmax=Tmax*tau
  if(gamma>0){
    tmp=lambda/(2*pi*C_0)*(sqrt(gamma)*atan(sqrt(gamma)*r/sqrt(2*D))/(2^1.5*D^1.5)+1/(2*D*r)-pi*sqrt(gamma)/(2^2.5*D^1.5))+1
  }else{ #This corresponds to the case U=0
        if(is.na(Tmax)){
                print("Tmax should have a value")
                stop()
        }
    tmp=lambda/(4*D*pi*C_0*r)*(1-erf(r/(sqrt(8*Tmax*D))))+1
  }
  return(tmp)
}

