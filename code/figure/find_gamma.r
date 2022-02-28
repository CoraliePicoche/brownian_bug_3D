#23/01/2021 CP: This script uses the output from  compute_gamma to compute the distance r(t) between pairs of particles with time, and then compute gamma as the slope of (1/d)<ln(r(t))> where <ln(r(t))> is the average obtained with all the 800 pairs of particles. We do that for 4 different values of Utot. This will be used for eq(2) for Fig.3

rm(list=ls())
graphics.off()

x=seq(-2*pi,2*pi,0.1)

#List of estimated coefficients
coef_estimated=c()

pdf("gamma_compute_lowest_adv.pdf")
par(mfrow=c(1,1))

f=read.table("../simulation/Spatial_distribution_201.txt",header=F,sep=";",dec='.')
colnames(f)=c("t","x","y","z","yfirst","type","first_parent")

f_param=read.table(paste("../simulation/param_201.txt",sep=""),header=F,sep="=",dec='.')
f_param[,2]=as.numeric(as.character(f_param[,2]))
tau=f_param[f_param[,1]=="tau",2]
Utot_tmp=f_param[f_param[,1]=="Utot",2]

print(paste("Utot",Utot_tmp))
timestep=unique(f$t)
unique_pair=unique(f$first_parent)

#To compute the linear coefficient, we need to only examine the separation as a function of time before steady state. We therefore only use a subset of timestep
	if (Utot_tmp=="0.5"){ #Stabilizes quickly: 15 time steps seems too be a good proxy
		subset=1:15
	}else if (Utot_tmp=="2.5"){ ##Stabilizes very quickly: 4 time steps and stable state
		subset=1:4
	}else{subset=1:length(timestep)}

#Compute distance between parent and children
	dist=rep(0,length(subset))
	dist[1]=log(10^(-7))*1/3 #The first distance between particle is  10-7
	for (p in 1:length(unique_pair)){ #A pair is defined by a unique id, pair, for particles identified as "P" or "C"
		pair=unique_pair[p] 
		f_pair=subset(f,first_parent==pair)
		for(t in 2:subset[length(subset)]){ #For each time step, we compute the mean distance between particles
			xP=f_pair$x[f_pair$t==timestep[t]&f_pair$type=="P"]
			yP=f_pair$y[f_pair$t==timestep[t]&f_pair$type=="P"]
			zP=f_pair$z[f_pair$t==timestep[t]&f_pair$type=="P"]
			xC=f_pair$x[f_pair$t==timestep[t]&f_pair$type=="C"]
			yC=f_pair$y[f_pair$t==timestep[t]&f_pair$type=="C"]
			zC=f_pair$z[f_pair$t==timestep[t]&f_pair$type=="C"]
			dist_tmp=sqrt((xP-xC)^2+(yP-yC)^2+(zP-zC)^2)
			if(dist_tmp==0.0){dist_tmp=10^(-9)} #If the distance is too small, we set an arbitrary low value
			if(is.nan(log(dist_tmp))){stop()}
			dist[t]=dist[t]+log(dist_tmp) #We sum all distances and will then divide by the number of particles
		}
	}
	dist[2:length(dist)]=(1/3)*dist[2:length(dist)]/length(unique_pair)

	modif_t=timestep[subset]*tau
	l=lm(dist~(modif_t)) #Finally, we estimate the value of the slope
	slope=l$coefficients[2]
	coef_estimated=c(coef_estimated,slope)
	plot(modif_t,dist,main=eval(substitute(expression(paste("U",tau,"/3=",Utot_tmp,", ",gamma,"=",slope,sep="")),list(Utot_tmp=Utot_tmp,slope=format(slope,digits=3)))),xlab="t",ylab="log(r)")
	lines(modif_t,l$coefficients[1]+modif_t*l$coefficient[2])

dev.off()
