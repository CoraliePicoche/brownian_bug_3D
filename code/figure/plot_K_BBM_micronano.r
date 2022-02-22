rm(list=ls())
graphics.off()

library(spatstat)
source("utilitary_functions.r")
source("theoretical_functions.r")

colo=c("red","blue","grey")

nb_simu_list=c(0,1,2,3)

a_gamma_list=c(1231,0,1231,0)

a_main_list=rep(c("Advection","No advection"),2)

a_text=c("a","b","c","d")

pdf("K_micronano.pdf",width=8)
par(mfrow=c(2,2))


for(i in 1:length(nb_simu_list)){

nb_simu=nb_simu_list[i]
a_gamma=a_gamma_list[i]

f=read.table(paste("../simulation/lambda_K_",nb_simu,".txt",sep=""),sep=";",header=F,dec=".")
colnames(f)=c("r","sp1","sp2","pcf","dominance","lambda_K","K")
unique_sp=unique(f$sp1)
if (length(unique_sp)>length(colo)){
        colo=rainbow(length(unique_sp))
}

#Count species
f_count=read.table(paste("../simulation/nb_indiv_",nb_simu,".txt",sep=""),header=F,sep=";",dec='.')
colnames(f_count)=c("species","abundance")

#Param simu
f_param=read.table(paste("../simulation/param_",nb_simu,".txt",sep=""),header=F,sep="=",dec='.')
f_param[,2]=as.numeric(as.character(f_param[,2]))
a_lambda=f_param[f_param[,1]=="growth_rate",2]
a_Delta=f_param[f_param[,1]=="Delta",2]
a_tau=f_param[f_param[,1]=="tau",2]
a_r=f_param[f_param[,1]=="radius",2]*10^2 #Radius is converted to cm here
lim=a_r*10*2 

#Spatial distribution
tab_repart=read.table(paste("../simulation/Spatial_distribution_",nb_simu,".txt",sep=""), header=F, sep=";")
colnames(tab_repart)=c("t","x","y","z","yfirst","ancestor","species")

Concentration=f_count$abundance/f_param[f_param[,1]=="volume",2]

if(i<=2){
	par(mar=c(3,6,4,2))
	xl=""
	ml=a_main_list[i]
}else{
	par(mar=c(4,6,3,2))
	xl="r"
	ml=""
}

if(i%%2==1){
	yl="K"
}else{
	yl=""
}

plot(0.1,0.1,log="xy",xlab=xl,ylab=yl,axes=F,xlim=range(f$r[f$r>0]),ylim=range(f$K[f$K>0]),cex.lab=1.5,t="n",main=ml)
axis(1, at=log10Tck('x','major'), tcl= 0.5,cex.axis=1.5) # bottom
axis(1, at=log10Tck('x','minor'), tcl= 0.1, labels=NA) # bottom
axis(2, at=log10Tck('y','major'), tcl= 0.5,cex.axis=1.5) # bottom
axis(2, at=log10Tck('y','minor'), tcl= 0.1, labels=NA) # bottom
box()
abline(v=lim,col="black",lty=4,lwd=2)
mtext(a_text[i],side=3,line=1.5,font=2,at=8*10^(-5))

if(i==1){
	mtext("Micro", outer = TRUE, cex = 1.5,line=-1.5,side=2,at=0.74)
	legend("topleft",c("Simu intra","Simu inter","Theory intra","Theory inter"),lty=c(1,NA,2,3),pch=c(NA,1,NA,NA),col=c("red","red","black","black"),bty="n",lwd=2)
}else if(i==3){
	mtext("Nano", outer = TRUE, cex = 1.5,line=-1.5,side=2,at=0.26)
}

for (sp in unique(f_count$species)){
	f_plot=subset(f,sp1==sp&sp2==sp)

#Ripley's function
	#Theory
	th_bbm=BBM_cdf(f_plot$r,a_gamma,a_Delta,Concentration[sp+1],a_lambda,a_tau,t=1000)
	th_poisson=4/3*pi*(f_plot$r)^3

	if(i==0){
	}

	for (spbis in unique(f_count$species)){
		if(spbis!=sp){
	        	f_plot_bis=subset(f,sp1==sp&sp2==spbis)
			points(f_plot_bis$r,f_plot_bis$K,col=colo[spbis],lwd=2)
		}
	}

	lines(f_plot$r,f_plot$K,col=colo[sp+1],lwd=2)
	lines(f_plot$r,th_bbm,lty=2,col="black",lwd=2)
	lines(f_plot$r,th_poisson,lty=3,col="black",lwd=2)
}
}
dev.off()
