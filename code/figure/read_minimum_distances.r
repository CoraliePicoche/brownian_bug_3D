rm(list=ls())
graphics.off()

nb_simu_list=c(152,153)
amain_list=c("Advection","No advection")
colo=c("red","blue","black")

pdf("distrib_distance_nano.pdf",width=7.5,height=10)
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
abline(h=0.554/(intensity^(1/3)),col="grey",lty=2,lwd=3)

for(a_sp in unique_sp){
	
	tab_1=subset(tab, sp==a_sp)
	dist_mono=tab_1[,3+a_sp]

	print(paste("For sp",a_sp))
	print(paste("MONO: Mean=",mean(dist_mono)," and min=",min(dist_mono),sep=""))
	points(a_sp,mean(dist_mono),pch=16,cex=1.5,col=colo[a_sp+1])
	points(a_sp,min(dist_mono),pch=17,col=colo[a_sp+1],cex=1.5)
	
	osp=setdiff(unique_sp,a_sp)
	for(b_sp in 1:length(osp)){
		dist_inter1=tab_1[,3+osp[b_sp]]
		points(a_sp+0.25*(b_sp/length(osp)),mean(dist_inter1),pch=1,cex=1.5,col=colo[osp[b_sp]+1])
		points(a_sp+0.25*(b_sp/length(osp)),min(dist_inter1),pch=2,col=colo[osp[b_sp]+1],cex=1.5)
		print(paste("CROSS:",osp[b_sp]," Mean=",mean(dist_inter1)," and min=",min(dist_inter1),sep=""))
	
#		lines(c(a_sp+1.15,a_sp+1.5),rep(0.554/(((realization[a_sp+1]+realization[osp[b_sp]+1])*intensity)^(1/3)),2),col="red",lwd=3,lty=1)

	}

	lines(c(a_sp,a_sp+0.25),rep(0.554/((realization[a_sp+1]*intensity)^(1/3)),2),col="green",lwd=3,lty=1)
} #end loop on unique_sp
if(i==1){
legend("topleft",c("Mean theory","Min theory","Monospecific","Mean Interspecific1","Min Interspecific2"),lty=c(2,1,NA,NA,NA),pch=c(NA,NA,16,1,2),col=c("grey","green","black","black","blue"),ncol=2,lwd=2)
mtext("Nano", outer = TRUE, cex = 1.5,line=-1.5)
} 
}

dev.off()
