rm(list=ls())
graphics.off()

source("utilitary_functions.r")
source("theoretical_functions.r")

nb_simu_tot=c(21,23)
nb_simu_na=list(c(25,27),26)
palette=c("cyan","blue","darkblue","darkcyan")

s=2 #We choose only one species

pdf("shift_of_dominance_noadv_zoom.pdf",width=13)
par(mfcol=c(1,2))


for(n in 1:length(nb_simu_tot)){
	    nb_simu=nb_simu_tot[n]

if(n==1){
	yl="dominance"
 	ml="Micro"
}else{
      	yl=""
       	ml="Nano"
}


#Let's add the dominance index for Thomas distribution
f_tot=read.table(paste("../simulation/lambda_K_",nb_simu,"_init.txt",sep=""),sep=";",header=F,dec=".")
colnames(f_tot)=c("r","sp1","sp2","pcf","dominance","lambda_K","K")
unique_sp=unique(f_tot$sp1)
u_r=unique(f_tot$r)

f_param=read.table(paste("../simulation/param_",nb_simu,".txt",sep=""),sep="=",header=F,dec=".")
f_param[,2]=as.numeric(as.character(f_param[,2]))
colnames(f_param)=c("name","value")
volume=f_param[f_param$name=="volume","value"]
N_parent=f_param[f_param$name=="N_parent","value"]
a_sigma=f_param[f_param$name=="sigma","value"]

concentrations=c(f_param[f_param$name=="init_size 0","value"],f_param[f_param$name=="init_size 1","value"],f_param[f_param$name=="init_size 2","value"])/volume

pow_th=seq(-4,0,length.out=50)
th_r=10^pow_th

th_poisson=4/3*pi*th_r^3
th_thomas=thomas_cdf(th_r,a_sigma,N_parent/volume)


plot(0,0,t="n",xlim=c(10^(-2),1),ylim=c(0.3,1),xlab="r",ylab=yl,log="x",main=ml,axes=F,cex.lab=1.5,cex.main=1.25)
axis(1, at=log10Tck('x','major'), tcl= 0.5,cex.axis=1.5) # bottom
axis(1, at=log10Tck('x','minor'), tcl= 0.1, labels=NA) # bottom
axis(2,tcl=0.5,cex.axis=1.5) # left
box()

th_dominance_thomas=concentrations[s]*(th_thomas)/(sum(concentrations)*th_poisson+concentrations[s]*(th_thomas-th_poisson))
f_tmp=subset(f_tot,sp1==unique_sp[s] & sp2==unique_sp[s])
id_b=floor(seq(1,length(f_tmp$r),length.out=50))
lines(f_tmp$r[id_b],f_tmp$dominance[id_b],col=palette[1],lty=1,lwd=3)
lines(th_r,th_dominance_thomas,col="black",lty=3,lwd=2)


#Then check after simulations
###WITH SHORT TIME

f_tot=read.table(paste("../simulation/lambda_K_",nb_simu,".txt",sep=""),sep=";",header=F,dec=".")
colnames(f_tot)=c("r","sp1","sp2","pcf","dominance","lambda_K","K")
unique_sp=unique(f_tot$sp1)
u_r=unique(f_tot$r)

f_count=read.table(paste("../simulation/nb_indiv_",nb_simu,".txt",sep=""),sep=";",header=F,dec=".")
colnames(f_count)=c("species","abundance")
concentrations=f_count$abundance/volume

Utot=f_param[f_param$name=="Utot","value"]
if(Utot==0.5){
	gamma=1231
}else if (Utot==0.0){
	gamma=0
}
a_Delta=f_param[f_param$name=="Delta","value"]
a_lambda=f_param[f_param$name=="growth_rate","value"]
a_tau=f_param[f_param$name=="tau","value"]
a_tmax=f_param[f_param$name=="tmax","value"]

th_bbm=BBM_cdf(th_r,gamma,a_Delta,concentrations[s],a_lambda,a_tau,t=a_tmax)
th_dominance=concentrations[s]*(th_bbm+th_thomas-th_poisson)/(sum(concentrations)*th_poisson+concentrations[s]*(th_bbm+th_thomas-2*th_poisson))
f_tmp=subset(f_tot,sp1==unique_sp[s] & sp2==unique_sp[s])
lines(f_tmp$r,f_tmp$dominance,col=palette[2])
points(th_r,th_dominance,col="black",pch=1)

for(m in 1:length(nb_simu_na[[n]])){
	print(nb_simu_na[[n]])
	print(nb_simu_na[[n]][m])
#Then check after simulations
##WITH LONG_TIME
f_tot=read.table(paste("../simulation/lambda_K_",nb_simu_na[[n]][m],".txt",sep=""),sep=";",header=F,dec=".")
colnames(f_tot)=c("r","sp1","sp2","pcf","dominance","lambda_K","K")
unique_sp=unique(f_tot$sp1)

u_r=unique(f_tot$r)

f_count=read.table(paste("../simulation/nb_indiv_",nb_simu_na[[n]][m],".txt",sep=""),sep=";",header=F,dec=".")
colnames(f_count)=c("species","abundance")
concentrations=f_count$abundance/volume
print(concentrations)

f_param=read.table(paste("../simulation/param_",nb_simu_na[[n]][m],".txt",sep=""),sep="=",header=F,dec=".")
f_param[,2]=as.numeric(as.character(f_param[,2]))
colnames(f_param)=c("name","value")

Utot=f_param[f_param$name=="Utot","value"]
if(Utot==0.5){
        gamma=1231
}else if (Utot==0.0){
        gamma=0
}
a_Delta=f_param[f_param$name=="Delta","value"]
a_lambda=f_param[f_param$name=="growth_rate","value"]
a_tau=f_param[f_param$name=="tau","value"]
a_tmax=f_param[f_param$name=="tmax","value"]

th_bbm=BBM_cdf(th_r,gamma,a_Delta,concentrations[s],a_lambda,a_tau,t=a_tmax)
th_dominance=concentrations[s]*(th_bbm+th_thomas-th_poisson)/(sum(concentrations)*th_poisson+concentrations[s]*(th_bbm+th_thomas-2*th_poisson))
f_tmp=subset(f_tot,sp1==unique_sp[s] & sp2==unique_sp[s])
lines(f_tmp$r,f_tmp$dominance,col=palette[2+m],lwd=2)
points(th_r,th_dominance,col="black",pch=c(0,2)[m])
}

if(n==1){
 	legend("bottomleft",c("t=0",expression('t='~10^3),expression('t='~10^5),expression('t='~10^6)),col=palette,bty="n",lwd=2,lty=1,pch=NA)
        mtext("a",side=3,line=1.5,font=2,at=1*10^(-2),cex=1.25)
}else{
        legend("bottomleft",c("Init Thomas",expression("Theory t="~10^3),expression("Theory t="~10^5),expression("Theory t="~10^6)),col="black",pch=c(NA,1,0,2),lty=c(3,NA,NA,NA),bty="n",lwd=2)
        mtext("b",side=3,line=1.5,font=2,at=1*10^(-2),cex=1.25)
}

}
dev.off()

