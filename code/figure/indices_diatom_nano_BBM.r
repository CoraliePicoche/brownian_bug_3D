rm(list=ls())
graphics.off()

library(pracma)

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

pdf("indices_diatom_nano_3sp_no_advection.pdf",width=10)
par(mfrow=c(2,3),mar=c(4,5,1,1))

#Using the delta=10^-5
#Diatom
#sim_diatom=list(29:32) #3 species
#sim_diatom=list(41:44) #10 species, skewed abundance
sim_diatom=list(60:63) #3 species, NO advection
#lim_min_diatom=5*10^(-3)
#lim_max_diatom=2.4*10^(-1) ##Phycosphere
lim_max_diatom=25*10^(-4)*10

#Nano
#sim_nano=list(37:40) #3 species
#sim_nano=list(45:48) #10 species, skewed abundance
sim_nano=list(64:67) #3 species, NO advection
#lim_min_nano=3*10^(-4) 
#lim_max_nano=10^(-3) ##Phycosphere
lim_max_nano=1.5*10^(-4)*10

tot_sim=list(sim_diatom,sim_nano)
#tot_min=list(lim_min_diatom,lim_min_nano)
tot_max=list(lim_max_diatom,lim_max_nano)

for(type in 1:length(tot_sim)){
	nb_simu_tot=tot_sim[type][[1]][[1]]

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
	        if (length(unique_sp)>length(colo)){
        	        colo=rainbow(length(unique_sp))
	        }
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

	if(type==2){
		xl="r"
	}else{
		xl=""
	}

        plot(0.1,0.1,t="n",log="xy",xlab=xl,ylab=expression(g[ii]~(r)),axes=F,xlim=range(f_tot$r[f_tot$r>0]),cex.lab=1.5,ylim=range(f_tot$pcf[f_tot$pcf>0],na.rm=T))
	#abline(v=c(tot_min[[type]][[1]],tot_max[[type]][[1]]),lty=3)
	#abline(v=c(tot_max[[type]][[1]]),lty=3)
	for (s1 in 1:length(unique_sp)){
        	f_plot=subset(f_tot,sp1==unique_sp[s1]&sp2==unique_sp[s1])
	        axis(1, at=log10Tck('x','major'), tcl= 0.5,cex.axis=1.5) # bottom
        	axis(1, at=log10Tck('x','minor'), tcl= 0.1, labels=NA) # bottom
	        axis(2, at=log10Tck('y','major'), tcl= 0.5,cex.axis=1.5) # bottom
        	axis(2, at=log10Tck('y','minor'), tcl= 0.1, labels=NA) # bottom
	        box()

                volume=unique(f_param_tot$value[f_param_tot$name=="volume"])
                lambda=unique(f_param_tot$value[f_param_tot$name=="growth_rate"])
                Delta=unique(f_param_tot$value[f_param_tot$name=="Delta"])
                T_max=unique(f_param_tot$value[f_param_tot$name=="tmax"])
                C_0=unique(f_count_tot$abundance[f_count_tot$species==unique_sp[s1]]/volume)
                U_tot=unique(f_param_tot$value[f_param_tot$name=="Utot"])

		if(U_tot==0){
			gamma=0
		}else if(U_tot==0.5){
			gamma=0.5
		}else{
			stop("you don't have the value of gamma")
		}

		s2=s1
                f_plot=subset(f_tot,sp1==unique_sp[s1]&sp2==unique_sp[s2])
        	points(f_plot$r,f_plot$pcf,lty=1,col=colo[s2])
                
		th_bbm=G_theoretical(gamma,f_plot$r,lambda,Delta,C_0,Tmax=T_max*0.0002)
                lines(f_plot$r,th_bbm)

		if(s1==3){
	       	 	plot(0.1,0.1,t="n",log="xy",xlab=xl,ylab=expression(g[2~j]~(r)),axes=F,xlim=range(f_tot$r[f_tot$r>0]),cex.lab=1.5,ylim=range(f_tot$pcf[f_tot$pcf>0],na.rm=T))
			axis(1, at=log10Tck('x','major'), tcl= 0.5,cex.axis=1.5) # bottom
                	axis(1, at=log10Tck('x','minor'), tcl= 0.1, labels=NA) # bottom
	                axis(2, at=log10Tck('y','major'), tcl= 0.5,cex.axis=1.5) # bottom
        	        axis(2, at=log10Tck('y','minor'), tcl= 0.1, labels=NA) # bottom
                	box()
			legend("topleft",paste("Sp=",unique_sp),col=colo[unique_sp],lty=1,bty="n")
		#abline(v=c(tot_min[[type]][[1]],tot_max[[type]][[1]]),lty=3)
		#abline(v=c(tot_max[[type]][[1]]),lty=3)

        		for(s2 in 1:length(unique_sp)){
                		f_plot=subset(f_tot,sp1==unique_sp[s1]&sp2==unique_sp[s2])
				if(s2!=s1){
			                points(f_plot$r,f_plot$pcf,lty=1,col=colo[s2],lwd=2)
				}
		        } #end loop on sp s2
	
		}#end test on s1
	} #end loop on species
	plot(1,1,t="n",xlim=range(f_tot$r[f$r>0]),ylim=range(f_tot$dominance,na.rm=T),xlab=xl,ylab="dominance",log="x",cex.lab=1.5,cex.axis=1.5,axes=F)
	axis(1, at=log10Tck('x','major'), tcl= 0.5,cex.axis=1.5) # bottom
	axis(1, at=log10Tck('x','minor'), tcl= 0.1, labels=NA) # bottom
	axis(2,tcl=0.5,cex.axis=1.5) # left
	box()

	for (s1 in 1:length(unique_sp)){
        	f_plot=subset(f_tot,sp1==unique_sp[s1]&sp2==unique_sp[s1])
	       	points(f_plot$r,f_plot$dominance,col=colo[s1])
	}

	#abline(v=c(tot_min[[type]][[1]],tot_max[[type]][[1]]),lty=3)
	#abline(v=tot_max[[type]][[1]],lty=3)
} #end loop on organism 

dev.off()
