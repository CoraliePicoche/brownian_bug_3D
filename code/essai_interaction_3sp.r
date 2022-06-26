rm(list=ls())
graphics.off()

set.seed(42)

library(deSolve)

state=c(N1=1,N2=1,N3=1)

nb_species=c()

alpha_ij=0.06
noise=0.001
#I either fix all parameters to be equal, or add some noise. I observe coexistence in the absence of noise, i.e. all interspecific interactions are the same

parameters=c(r=1,alpha_ii=1,
	     alpha_12=alpha_ij+rnorm(1,0,noise),
	     alpha_13=alpha_ij+rnorm(1,0,noise),
	     alpha_21=alpha_ij+rnorm(1,0,noise),
	     alpha_23=alpha_ij+rnorm(1,0,noise),
	     alpha_31=alpha_ij+rnorm(1,0,noise),
	     alpha_32=alpha_ij+rnorm(1,0,noise))

times=seq(1,10000,by=0.01)

func_LV=function(times,state,parameters)
	with(as.list(c(state, parameters)), {
		dN1 <- N1*(r-alpha_ii*N1-alpha_12*(N2)-alpha_13*N3)
		dN2 <- N2*(r-alpha_ii*N2-alpha_21*(N1)-alpha_23*N3)
		dN3 <- N3*(r-alpha_ii*N3-alpha_31*(N1)-alpha_32*N2)
		list(c(dN1, dN2, dN3))
})

out <- ode(y = state, times = times, func = func_LV, parms = parameters)
plot(out)
