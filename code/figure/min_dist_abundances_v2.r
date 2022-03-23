rm(list=ls())
graphics.off()

nb_simu_list=c(10,12)

pdf("dist_abundances_10sp_v2.pdf")

lims=c(0.01,0.8)
plot(0,0,t="n",xlim=c(500,5e6),ylim=lims,log="x",xlab="C/L",ylab="Mean dist. to nearest neighbour",xaxt="n")
axis(1)

for(i in 1:length(nb_simu_list)){
	nb_simu=nb_simu_list[i]

	if(i==1){colo="black"}
	else{colo="grey"}

tab=read.table(paste("../simulation/min_dist_",nb_simu,".txt",sep=""),header=F,sep=";")
colnames(tab)=c("p","sp","sp1","sp2","sp3")

f_param=read.table(paste("../simulation/param_",nb_simu,".txt",sep=""),header=F,sep="=")
f_param[,2]=as.numeric(as.character(f_param[,2]))
volume=f_param[f_param[,1]=="volume",2]*10^(-3)

nb_indiv=read.table(paste("../simulation/nb_indiv_",nb_simu,".txt",sep=""),header=F,sep=";")
realization=nb_indiv[,2]

unique_sp=unique(tab$sp)


for(a_sp in unique_sp){
	
	tab_1=subset(tab, sp==a_sp)
	dist_mono=tab_1[,3+a_sp]

	print(paste("For sp",a_sp))
	print(paste("MONO: Mean=",mean(dist_mono)," and min=",min(dist_mono),sep=""))
	points(realization[a_sp+1]/volume,mean(dist_mono),pch=16,cex=1.,col=colo)
	
	osp=setdiff(unique_sp,a_sp)
	dist=c()
	for(b_sp in 1:length(osp)){
		dist=c(dist,tab_1[,3+osp[b_sp]])
	}
	points(realization[a_sp+1]/volume,mean(dist),pch=0,col=colo)

} #end loop on unique_sp
}	

legend("topright",c("Intraspecific","Interspecific","Micro","Nano"),pch=c(16,0,16,16),col=c("black","black","black","grey"),lwd=2,bty="n",lty=NA)
dev.off()
