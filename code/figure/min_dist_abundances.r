rm(list=ls())
graphics.off()

nb_simu_list=c(10,11,12,13)
amain_list=c("Advection","No advection")
colo=c("red","blue","black")
colo=rainbow(10)

pdf("dist_abundances_10sp.pdf",width=7.5,height=10)
par(mfrow=c(2,2))

for(i in 1:length(nb_simu_list)){
	nb_simu=nb_simu_list[i]

tab=read.table(paste("../simulation/min_dist_",nb_simu,".txt",sep=""),header=F,sep=";")
colnames(tab)=c("p","sp","sp1","sp2","sp3")

f_param=read.table(paste("../simulation/param_",nb_simu,".txt",sep=""),header=F,sep="=")
f_param[,2]=as.numeric(as.character(f_param[,2]))
volume=f_param[f_param[,1]=="volume",2]*10^(-3)

nb_indiv=read.table(paste("../simulation/nb_indiv_",nb_simu,".txt",sep=""),header=F,sep=";")
realization=nb_indiv[,2]

unique_sp=unique(tab$sp)

#plot(0,0,t="n",xlim=range(realization),ylim=range(tab[,3:ncol(tab)],na.rm=T),log="xy",xlab="Species",ylab="Dist. to nearest neighbour",xaxt="n",main=amain_list[i])
if(i<3){
	lims=c(0.1,0.8)
}else{
	lims=c(0.01,0.2)
}
plot(0,0,t="n",xlim=range(realization)/(volume),ylim=lims,log="x",xlab="C/L",ylab="Mean dist. to nearest neighbour",xaxt="n",main=amain_list[i])
axis(1)

for(a_sp in unique_sp){
	
	tab_1=subset(tab, sp==a_sp)
	dist_mono=tab_1[,3+a_sp]

	print(paste("For sp",a_sp))
	print(paste("MONO: Mean=",mean(dist_mono)," and min=",min(dist_mono),sep=""))
#	boxplot(dist_mono,at=realization[a_sp+1],col=colo[a_sp+1],add=T,border=colo[a_sp+1],pt.col=colo[a_sp+1],range=0)
	points(realization[a_sp+1]/volume,mean(dist_mono),pch=16,cex=1.,col="black")
#	points(a_sp,min(dist_mono),pch=24,bg=colo[a_sp+1],cex=1.5,col="black")
	
	osp=setdiff(unique_sp,a_sp)
	dist=c()
#	list_dist_mean=rep(NA,length(osp))
#	list_dist_min=rep(NA,length(osp))
	for(b_sp in 1:length(osp)){
		dist=c(dist,tab_1[,3+osp[b_sp]])
#		list_dist_mean[b_sp]=mean(tab_1[,3+osp[b_sp]])
#		list_dist_min[b_sp]=min(tab_1[,3+osp[b_sp]])
#		print(paste("CROSS:",osp[b_sp]," Mean=",list_dist_mean[b_sp]," and min=",list_dist_min[b_sp],sep=""))
	}
	points(realization[a_sp+1]/volume,mean(dist),pch=0,col="black")
#	boxplot(list_dist_mean,at=a_sp+0.35,col=colo[a_sp+1],add=T,border=colo[a_sp+1],pt.col=colo[a_sp+1])
#	boxplot(list_dist_min,at=a_sp+0.35,col=colo[a_sp+1],add=T,border="black")
#	boxplot(dist,at=10+realization[a_sp+1],col=colo[a_sp+1],add=T,border="black",range=0)

#	lines(c(a_sp,a_sp+0.25),rep(0.554/((realization[a_sp+1]*intensity)^(1/3)),2),col="black",lwd=3,lty=1)
} #end loop on unique_sp
if(i==1){
legend("topright",c("Intraspecific","Interspecific"),pch=c(16,0),col=c("black"),lwd=2,bty="n",lty=NA)
mtext("Micro", outer = TRUE, cex = 1.5,line=-1.5)
}else if(i==3){
mtext("Nano", outer = TRUE, cex = 1.5,line=-32)
}	
}

dev.off()
