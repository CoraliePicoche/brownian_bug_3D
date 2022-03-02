rm(list=ls())
graphics.off()

nb_simu_list=c(30,31)
amain_list=c("Advection","No advection")
colo=c("red","blue","black")
colo=rainbow(10)

pdf("distrib_distance_micro_box_10sp.pdf",width=7.5,height=10)
par(mfrow=c(2,1))

for(i in 1:length(nb_simu_list)){
	nb_simu=nb_simu_list[i]

tab=read.table(paste("../simulation/min_dist_",nb_simu,".txt",sep=""),header=F,sep=";")
colnames(tab)=c("p","sp","sp1","sp2","sp3")

f_param=read.table(paste("../simulation/param_",nb_simu,".txt",sep=""),header=F,sep="=")
f_param[,2]=as.numeric(as.character(f_param[,2]))
intensity=f_param[f_param[,1]=="init_size 0",2]/f_param[f_param[,1]=="volume",2]

nb_indiv=read.table(paste("../simulation/nb_indiv_",nb_simu,".txt",sep=""),header=F,sep=";")
realization=nb_indiv[,2]

unique_sp=unique(tab$sp)

plot(0,0,t="n",xlim=c(0,max(unique_sp)+0.5),ylim=c(10^(-5),10^(0)),log="y",xlab="Species",ylab="Dist. to nearest neighbour",xaxt="n",main=amain_list[i])
axis(1,at=unique_sp)
abline(h=0.554/(intensity^(1/3)),col="grey",lty=1,lwd=3)

for(a_sp in unique_sp){
	
	tab_1=subset(tab, sp==a_sp)
	dist_mono=tab_1[,3+a_sp]

	print(paste("For sp",a_sp))
	print(paste("MONO: Mean=",mean(dist_mono)," and min=",min(dist_mono),sep=""))
	points(a_sp,mean(dist_mono),pch=21,cex=1.5,bg=colo[a_sp+1],col="grey")
	points(a_sp,min(dist_mono),pch=24,bg=colo[a_sp+1],cex=1.5,col="black")
	
	osp=setdiff(unique_sp,a_sp)
	list_dist_mean=rep(NA,length(osp))
	list_dist_min=rep(NA,length(osp))
	for(b_sp in 1:length(osp)){
		list_dist_mean[b_sp]=mean(tab_1[,3+osp[b_sp]])
		list_dist_min[b_sp]=min(tab_1[,3+osp[b_sp]])
		print(paste("CROSS:",osp[b_sp]," Mean=",list_dist_mean[b_sp]," and min=",list_dist_min[b_sp],sep=""))
	}
	boxplot(list_dist_mean,at=a_sp+0.35,col=colo[a_sp+1],add=T,border=colo[a_sp+1],pt.col=colo[a_sp+1])
	boxplot(list_dist_min,at=a_sp+0.35,col=colo[a_sp+1],add=T,border="black")

	lines(c(a_sp,a_sp+0.25),rep(0.554/((realization[a_sp+1]*intensity)^(1/3)),2),col="black",lwd=3,lty=1)
} #end loop on unique_sp
if(i==1){
legend("bottomright",c("Mean","Minimum","Monospecific","Interspecific"),lty=c(1,1,NA,NA,NA),pch=c(NA,NA,17,15),col=c("grey","black","red","red"),lwd=2,pt.cex=c(1,1,1,2),ncol=2)
mtext("Micro", outer = TRUE, cex = 1.5,line=-1.5)
} 
}

dev.off()
