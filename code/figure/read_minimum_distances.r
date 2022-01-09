rm(list=ls())
graphics.off()

nb_simu=84
tab=read.table(paste("../simulation/min_dist_",nb_simu,".txt",sep=""),header=F,sep=";")
colnames(tab)=c("p","sp","sp1","sp2","sp3")

f_param=read.table(paste("../simulation/param_",nb_simu,".txt",sep=""),header=F,sep="=")
f_param[,2]=as.numeric(as.character(f_param[,2]))
intensity=f_param[f_param[,1]=="init_size 0",2]/f_param[f_param[,1]=="volume",2]

nb_indiv=read.table(paste("../simulation/nb_indiv_",nb_simu,".txt",sep=""),header=F,sep=";")
realization=nb_indiv[,2]

unique_sp=unique(tab$sp)

pdf("distrib_distance.pdf")
plot(0,0,t="n",xlim=c(1,3.5),ylim=c(10^(-5),10^(-1)),log="y",xlab="Species",ylab="Min distance")
abline(h=0.554/(intensity^(1/3)),col="grey",lty=2,lwd=3)

for(a_sp in unique_sp){

	tab_1=subset(tab, sp==a_sp)
	dist_mono=tab_1[,3+a_sp]
	osp=setdiff(unique_sp,a_sp)
	dist_inter1=tab_1[,3+osp[1]]
	dist_inter2=tab_1[,3+osp[2]]

	print(paste("For sp",a_sp))
	print(paste("MONO: Mean=",mean(dist_mono)," and min=",min(dist_mono),sep=""))
	print(paste("CROSS1: Mean=",mean(dist_inter1)," and min=",min(dist_inter1),sep=""))
	print(paste("CROSS2: Mean=",mean(dist_inter2)," and min=",min(dist_inter2),sep=""))

	points(a_sp+1,mean(dist_mono),pch=16,cex=1.5)
	points(a_sp+1+0.05,mean(dist_inter1),pch=1,cex=1.5)
	points(a_sp+1+0.1,mean(dist_inter2),pch=2,cex=1.5)
	
	points(a_sp+0.25+1,min(dist_mono),pch=16,col="blue",cex=1.5)
	points(a_sp+1+0.25+0.05,min(dist_inter1),pch=1,col="blue",cex=1.5)
	points(a_sp+1+0.25+0.1,min(dist_inter2),pch=2,col="blue",cex=1.5)
	
	lines(c(a_sp+1.25,a_sp+1.35),rep(0.554/((realization[a_sp+1]*intensity)^(1/3)),2),col="blue",lwd=3,lty=4)
}

legend("topleft",c("Mean theory","Min theory","Monospecific","Mean Interspecific1","Min Interspecific2"),lty=c(2,4,NA,NA,NA),pch=c(NA,NA,16,1,2),col=c("grey","blue","black","black","blue"),ncol=3,lwd=2)

dev.off()
