rm(list=ls())
graphics.off()

source("utilitary_functions.r")
source("theoretical_functions.r")

nb_simu_tot=c(20,22)
palette=c("red","blue","grey")

pdf("shift_of_dominance.pdf",width=13)
par(mfrow=c(1,2))

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

id=seq(1,length(u_r),length.out=50)

th_poisson=4/3*pi*u_r[id]^3
th_thomas=thomas_cdf(u_r[id],a_sigma)/(N_parent/volume)
plot(0,0,t="n",xlim=range(u_r),ylim=c(0.3,1),xlab="r",ylab=yl,log="x",main=ml,axes=F)
                        axis(1, at=log10Tck('x','major'), tcl= 0.5,cex.axis=1.5) # bottom
                        axis(1, at=log10Tck('x','minor'), tcl= 0.1, labels=NA) # bottom
                        axis(2,tcl=0.5,cex.axis=1.5) # left
			box()

for (s in 1:length(unique_sp)){
        th_dominance=concentrations[s]*(th_poisson+th_thomas)/(sum(concentrations)*th_poisson+concentrations[s]*th_thomas)
        f_tmp=subset(f_tot,sp1==unique_sp[s] & sp2==unique_sp[s])
        id_b=seq(1,length(f_tmp$r),length.out=50)
        points(f_tmp$r[id_b],f_tmp$dominance[id_b],col=palette[s],pch=2,cex=1)
        points(u_r[id],th_dominance,col="black",pch=2,cex=1)
}




#Then check after simulations
###WITH ADVECTION

f_tot=read.table(paste("../simulation/lambda_K_",nb_simu,".txt",sep=""),sep=";",header=F,dec=".")
colnames(f_tot)=c("r","sp1","sp2","pcf","dominance","lambda_K","K")
unique_sp=unique(f_tot$sp1)
u_r=unique(f_tot$r)
th_poisson=4/3*pi*u_r^3

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


for (s in 1:length(unique_sp)){
	th_bbm=BBM_cdf(u_r,gamma,a_Delta,concentrations[s],a_lambda,a_tau,t=a_tmax)
        th_dominance=concentrations[s]*(th_bbm)/(sum(concentrations)*th_poisson+concentrations[s]*(th_bbm-th_poisson))
        f_tmp=subset(f_tot,sp1==unique_sp[s] & sp2==unique_sp[s])
        points(f_tmp$r,f_tmp$dominance,col=palette[s],pch=16)
        points(u_r,th_dominance,col="black",pch=1)
}

#Then check after simulations
##WITHOUT ADVECTION
f_tot=read.table(paste("../simulation/lambda_K_",nb_simu+1,".txt",sep=""),sep=";",header=F,dec=".")
colnames(f_tot)=c("r","sp1","sp2","pcf","dominance","lambda_K","K")
unique_sp=unique(f_tot$sp1)

u_r=unique(f_tot$r)

f_count=read.table(paste("../simulation/nb_indiv_",nb_simu+1,".txt",sep=""),sep=";",header=F,dec=".")
colnames(f_count)=c("species","abundance")
concentrations=f_count$abundance/volume

f_param=read.table(paste("../simulation/param_",nb_simu+1,".txt",sep=""),sep="=",header=F,dec=".")
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

for (s in 1:length(unique_sp)){
        th_bbm=BBM_cdf(u_r,gamma,a_Delta,concentrations[s],a_lambda,a_tau,t=a_tmax)
        th_dominance=concentrations[s]*(th_bbm)/(sum(concentrations)*th_poisson+concentrations[s]*(th_bbm-th_poisson))
        f_tmp=subset(f_tot,sp1==unique_sp[s] & sp2==unique_sp[s])
        lines(f_tmp$r,f_tmp$dominance,col=palette[s],lwd=3)
        lines(u_r,th_dominance,col="black",lwd=3,lty=2)
}

if(n==1){
 	legend("bottomleft",c(paste("Sp=",unique_sp),"Theory"),col=c(palette[unique_sp+1],"black"),pch=1,bty="n",lwd=2,lty=NA)
        mtext("a",side=3,line=1.5,font=2,at=8*10^(-5))
}else{
        legend("bottomleft",c(expression(U~tau~"/"~3~"="~0),expression(U~tau~"/"~3~"="~0.5),"Init Thomas"),col="black",pch=c(NA,1,2),lty=c(1,NA,NA),bty="n",lwd=2)
        mtext("b",side=3,line=1.5,font=2,at=8*10^(-5))
}



}

dev.off()

