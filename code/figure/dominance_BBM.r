rm(list=ls())
graphics.off()

############################### Utilitary functions

#Ticks for log scale from https://stackoverflow.com/questions/47890742/logarithmic-scale-plot-in-r
log10Tck <- function(side, type){
  lim <- switch(side,
                x = par('usr')[1:2],
                y = par('usr')[3:4],
                stop("side argument must be 'x' or 'y'"))
  at <- floor(lim[1]) : ceiling(lim[2])
  return(switch(type,
                minor = outer(1:9, 10^(min(at):max(at))),
                major = 10^at,
                stop("type argument must be 'major' or 'minor'")
  ))
}

#Theoretical value of G
G_theoretical=function(gamma,r,lambda,Delta,C_0,Tmax=NA){ #U=0 in the absence of advection
  tau=0.0002
  D=Delta^2/(2*tau)
  if(gamma>0){
    tmp=lambda/(4*pi*C_0)*(sqrt(gamma)*atan(sqrt(gamma)*r/sqrt(2*D))/(2^1.5*D^1.5)+1/(2*D*r)-pi*sqrt(gamma)/(2^2.5*D^1.5))+1
  }else{ #This corresponds to the case U=0
        if(is.na(Tmax)){
                print("Tmax should have a value")
                stop()
        }
    tmp=2*lambda/(8*D*pi*C_0*r)*(1-erf(r/(2^1.5*sqrt(Tmax*D))))+1
  }
  return(tmp)
}

#####################################

colo=c("red","blue","grey","orange","violet")

pdf("dominance_BBM_diatom.pdf")
par(mar=c(4,4,1,1))

sim=list(21:24,25:28,29:32,33:36) #Diatom
lim_min=5*10^(-3)
lim_max=2.5*10^(-1)
#sim=list(37:40) #Nano
#lim_min=3*10^(-4)
#lim_max=10^(-2)
s1=1
legend_m=c()

plot(0.1,0.1,t="n",log="xy",xlab="r",ylab="g(r)",axes=F,xlim=c(10^(-4),1),cex.lab=1.5,ylim=c(0.3,1))
axis(1, at=log10Tck('x','major'), tcl= 0.5,cex.axis=1.5) # bottom
axis(1, at=log10Tck('x','minor'), tcl= 0.1, labels=NA) # bottom
axis(2,tcl=0.5,cex.axis=1.5) # left
box()
f_param_tot=matrix(,0,2)
colnames(f_param_tot)=c("name","value") #We initialize this one before the loop to have *all* values of parameters so that we can check they are all the same across the simulations
f_count_tot=matrix(,0,2)
colnames(f_count_tot)=c("species","abundance") #Same reason as above

for (s in 1:length(sim)){

	nb_simu_tot=sim[s][[1]]

	f_tot=matrix(,0,5)
	colnames(f_tot)=c("r","sp1","sp2","pcf","dominance")

	for (nb_simu in nb_simu_tot){
		f=read.table(paste("../simulation/pcf_",nb_simu,".txt",sep=""),sep=";",header=F,dec=".")
		colnames(f)=c("r","sp1","sp2","pcf","dominance")
		f_tot=rbind(f_tot,f)
		unique_sp=unique(f_tot$sp1)

		#Count species
		f_count=read.table(paste("../simulation/nb_indiv_",nb_simu,".txt",sep=""),header=F,sep=";",dec='.')
		colnames(f_count)=c("species","abundance")
		f_count_tot=rbind(f_count_tot,f_count)

		#Param simu
		f_param=read.table(paste("../simulation/param_",nb_simu,".txt",sep=""),header=F,sep="=",dec='.')
		f_param[,2]=as.numeric(as.character(f_param[,2]))
		colnames(f_param)=c("name","value")
		f_param_tot=rbind(f_param_tot,f_param)
	} #end loop on nb_simu
	delta_val=unique(f_param[f_param$name=="delta","value"])
	print(delta_val)
        if(s==1){
                m=paste(expression(delta),"=", delta_val, " spatstat",sep="")
        }else{
                m=paste(expression(delta), "=",delta_val,sep="")
        }
	legend_m=c(legend_m,m)

	f_plot=subset(f_tot,sp1==unique_sp[s1]&sp2==unique_sp[s1])

	points(f_plot$r,f_plot$dominance,lty=1,col=colo[s])
	lines(f_plot$r,f_plot$dominance,lty=1,col=colo[s])
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

nb_indiv=f_count_tot$abundance[f_count_tot$species==unique_sp[s1]]
if(length(unique(nb_indiv))>1){
	stop("nb indiv were different")
}
nb_indiv=unique(nb_indiv)
conc=nb_indiv/volume

if(U_tot==0.5){
	gamma=0.5
}else{
	stop("You don't have the right gamma")
}

legend("bottomleft",legend_m,lty=1,col=colo,bty="n",lwd=2)

abline(v=c(lim_min,lim_max),lty=3)

dev.off()
