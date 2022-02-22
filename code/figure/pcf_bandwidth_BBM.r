rm(list=ls())
graphics.off()

source("utilitary_functions.r")
source("theoretical_functions.r")

colo=c("darkorchid","green","darkblue","orange")

gamma=1231

pdf("bandwidth_BBM.pdf")
par(mar=c(4,4,1,1))

sim=list(4,5,6,0) #Diatoms
s1=1

plot(0.1,0.1,t="n",log="xy",xlab="r",ylab="g(r)",axes=F,xlim=c(10^(-4),1),cex.lab=1.5,ylim=c(1,2*10^7))
axis(1, at=log10Tck('x','major'), tcl= 0.5,cex.axis=1.5) # bottom
axis(1, at=log10Tck('x','minor'), tcl= 0.1, labels=NA) # bottom
axis(2, at=log10Tck('y','major'), tcl= 0.5,cex.axis=1.5) # bottom
axis(2, at=log10Tck('y','minor'), tcl= 0.1, labels=NA) # bottom
box()

f_param_tot=matrix(,0,2)
colnames(f_param_tot)=c("name","value") #We initialize this one before the loop to have *all* values of parameters so that we can check they are all the same across the simulations
f_count_tot=matrix(,0,2)
colnames(f_count_tot)=c("species","abundance") #Same reason as above

legend_m=c()

for (nb_simu in 1:length(sim)){
	f=read.table(paste("../simulation/lambda_K_",sim[nb_simu],".txt",sep=""),sep=";",header=F,dec=".")
	colnames(f)=c("r","sp1","sp2","pcf","dominance","lambda_K","K")
	unique_sp=unique(f$sp1)

	#Count species
	f_count=read.table(paste("../simulation/nb_indiv_",sim[nb_simu],".txt",sep=""),header=F,sep=";",dec='.')
	colnames(f_count)=c("species","abundance")
	f_count_tot=rbind(f_count_tot,f_count)

	#Param simu
	f_param=read.table(paste("../simulation/param_",sim[nb_simu],".txt",sep=""),header=F,sep="=",dec='.')
	f_param[,2]=as.numeric(as.character(f_param[,2]))
	colnames(f_param)=c("name","value")
	f_param_tot=rbind(f_param_tot,f_param)

delta_val=unique(f_param[f_param$name=="delta","value"])
legend_m=c(legend_m,eval(substitute(expression(paste(delta, "=",delta_val,sep="")),list(delta_val=delta_val))))

f_plot=subset(f,sp1==unique_sp[s1]&sp2==unique_sp[s1])

lines(f_plot$r,f_plot$pcf,lty=1,col=colo[nb_simu],t="o")
print(colo[nb_simu])
} #end loop on series of simulations


U_tot=f_param_tot[f_param_tot$name=="Utot","value"]	
if(length(unique(U_tot))>1){
	stop("sigma were different")
}
U_tot=unique(U_tot)

volume=f_param_tot[f_param_tot$name=="volume","value"]
if(length(unique(volume))>1){
	stop("volume were different")
}
volume=unique(volume)

Delta=f_param_tot[f_param_tot$name=="Delta","value"]
if(length(unique(Delta))>1){
	stop("Delta were different")
}
Delta=unique(Delta)
	
lambda=f_param_tot[f_param_tot$name=="growth_rate","value"]
if(length(unique(lambda))>1){
	stop("lambda were different")
}
lambda=unique(lambda)

tau=f_param_tot[f_param_tot$name=="tau","value"]
if(length(unique(tau))>1){
	stop("tau were different")
}
lambda=unique(lambda)

nb_indiv=f_count_tot$abundance[f_count_tot$species==unique_sp[s1]]
if(length(unique(nb_indiv))>1){
	stop("nb indiv were different")
}
nb_indiv=unique(nb_indiv)
conc=nb_indiv/volume

th_bbm=G_theoretical(gamma,unique(f$r),lambda,Delta,conc,tau,Tmax=NA)
lines(unique(f$r),th_bbm,lty=2,lwd=2)
legend("topright",legend_m,lty=1,col=colo,bty="n",lwd=2)

dev.off()
