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

#####################################

colo=c("red","blue","grey")

pdf("dominance_diatom_nano_compare_advection.pdf",width=10)
par(mfrow=c(1,2))

#Using the delta=10^-5
#Diatom
sim_diatom=list(29:32,60:63) #Advection, no advection
lim_max_diatom=25*10^(-4)*10

#Nano
sim_nano=list(37:40,64:67) #10 species
lim_max_nano=1.5*10^(-4)*10

tot_sim=list(sim_diatom,sim_nano)
tot_max=list(lim_max_diatom,lim_max_nano)

for(orga in 1:length(tot_sim)){ #Organism: diatom or nano
	for(adv in 1:length(tot_sim[orga][[1]])){ #Advection or no advection
		nb_simu_tot=tot_sim[orga][[1]][[adv]]
	
		f_tot=matrix(,0,5)
		colnames(f_tot)=c("r","sp1","sp2","pcf","dominance")
		f_count_tot=matrix(,0,2)
		colnames(f_count_tot)=c("species","abundance")
		f_param_tot=matrix(,0,2)
		colnames(f_param_tot)=c("name","value")

		for (nb_simu in nb_simu_tot){
        		f=read.table(paste("../simulation/pcf_",nb_simu,".txt",sep=""),sep=";",header=F,dec=".")
		        colnames(f)=c("r","sp1","sp2","pcf","dominance")
        		unique_sp=unique(f$sp1)
        		f_tot=rbind(f_tot,f)

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

		if(orga==1){
			yl="dominance"
			ml="Diatom"
		}else{
			yl=""
			ml="Nano"
		}

		if(adv==1){
			plot(1,1,t="n",xlim=range(f_tot$r[f$r>0]),ylim=range(f_tot$dominance,na.rm=T),xlab="r",ylab=yl,log="x",cex.lab=1.5,cex.axis=1.5,axes=F,main=ml)
			axis(1, at=log10Tck('x','major'), tcl= 0.5,cex.axis=1.5) # bottom
			axis(1, at=log10Tck('x','minor'), tcl= 0.1, labels=NA) # bottom
			axis(2,tcl=0.5,cex.axis=1.5) # left
			box()
		}

		if(orga==1){
			legend("bottomleft",paste("Sp=",unique_sp),col=colo[unique_sp+1],pch=1,bty="n",lwd=2,lty=NA)
		}else{
			legend("topright",c(expression(U~tau~"/"~2~"="~0),expression(U~tau~"/"~2~"="~0.5)),col="black",pch=c(NA,1),lty=c(1,NA),bty="n",lwd=2)
		}

		for (s1 in 1:length(unique_sp)){
        		f_plot=subset(f_tot,sp1==unique_sp[s1]&sp2==unique_sp[s1])
			if(adv==2){ #No advection
				if(orga==1){
		       			lines(f_plot$r,f_plot$dominance,col=colo[s1])
				}else{
		       			lines(f_plot$r,f_plot$dominance,col=colo[s1],pch=4)
				}
			}else{
				if(orga==1){
			       		points(f_plot$r,f_plot$dominance,col=colo[s1])
				}else{
		       			points(f_plot$r,f_plot$dominance,col=colo[s1],pch=4)
				}
			}
		}

	abline(v=c(tot_max[[orga]][[1]]),lty=2,col=c("green","darkblue")[orga],lwd=3)

	} #end loop on advection
} #end loop on organism

dev.off()
