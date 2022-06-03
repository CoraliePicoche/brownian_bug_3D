rm(list=ls())
graphics.off()

set.seed(42)

library(deSolve)

state=c(N1=1,N2=1)

nb_species=c()

alpha_ij=2
noise=0.0000
#I either fix all parameters to be equal, or add some noise. I observe coexistence in the absence of noise, i.e. all interspecific interactions are the same

parameters=c(r1=2,r2=2.1,alpha_ii=1,
	     alpha_12=alpha_ij+rnorm(1,0,noise),
	     alpha_21=alpha_ij+rnorm(1,0,noise))

print(parameters["alpha_ii"]^2)
print(parameters["alpha_12"]*parameters["alpha_21"])

times=seq(1,500,by=0.01)

func_LV=function(times,state,parameters)
	with(as.list(c(state, parameters)), {
		dN1 <- N1*(r1-alpha_ii*N1-alpha_12*(N2))#+N3))
		dN2 <- N2*(r2-alpha_ii*N2-alpha_21*(N1))#+N3))
		list(c(dN1, dN2))
})

out <- ode(y = state, times = times, func = func_LV, parms = parameters)
plot(out)
