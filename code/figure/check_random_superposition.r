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

thomas_cdf=function(r){
	sigma=0.01
        a=2/(sigma*sqrt(pi))*(erf(r/(2*sigma))-r/(sigma*sqrt(pi))*exp(-(r^2/(4*sigma^2))))
        return(a)
}


#####################################

colo=c("red","blue","grey")

pdf("K_Thomas.pdf",width=15)
par(mfrow=c(1,3),mar=c(4,4.5,3,1))

sim=list(101)


for (species in 1:3){
	s1=species
	ss=s1
	if(s1==1){
		yl="K(r)"
	}else{
		yl=""
	}
plot(0.1,0.1,t="n",log="x",xlab="r",ylab=yl,axes=F,xlim=c(10^(-4),0.1),cex.lab=1.5,ylim=c(10^(-6),(4/3)*pi*0.001),main=paste("Species=",s1,sep=""))
#plot(0.1,0.1,t="n",log="xy",xlab="r",ylab="K(r)",axes=F,xlim=c(10^(-4),1),cex.lab=1.5,ylim=c(10^(-7),10^(-2)))
axis(1, at=log10Tck('x','major'), tcl= 0.5,cex.axis=1.5) # bottom
axis(1, at=log10Tck('x','minor'), tcl= 0.1, labels=NA) # bottom
#axis(2, at=log10Tck('y','major'), tcl= 0.5,cex.axis=1.5) # bottom
#axis(2, at=log10Tck('y','minor'), tcl= 0.1, labels=NA) # bottom
axis(2,cex.axis=1.5)
box()
for (s in 1:length(sim)){

	nb_simu_tot=sim[s][[1]]

	f_tot=matrix(,0,5)
	colnames(f_tot)=c("r","sp1","sp2","pcf","dominance")
	f_count_tot=matrix(,0,2)
	colnames(f_count_tot)=c("species","abundance")
	f_param_tot=matrix(,0,2)
	colnames(f_param_tot)=c("name","value")

	for (nb_simu in nb_simu_tot){
		f=read.table(paste("../simulation/lambda_K_",nb_simu,".txt",sep=""),sep=";",header=F,dec=".")
		colnames(f)=c("r","sp1","sp2","pcf","dominance","lambda_K")
		f_tot=rbind(f_tot,f)
		unique_sp=unique(f_tot$sp1)

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
	sigma=f_param_tot[f_param_tot$name=="sigma","value"]
	if(length(unique(sigma))>1){
		stop("sigma were different")
	}
	sigma=unique(sigma)

	N_parent=f_param_tot[f_param_tot$name=="N_parent","value"]
	if(length(unique(N_parent))>1){
		stop("N_parent were different")
	}
	N_parent=unique(N_parent)

	volume=f_param_tot[f_param_tot$name=="volume","value"]
	if(length(unique(volume))>1){
		stop("volume were different")
	}
	volume=unique(volume)
	th_poisson=(4/3)*pi*unique(f_tot$r)^3

	f_plot=subset(f_tot,sp1==unique_sp[s1]&sp2==unique_sp[s1])

	points(f_plot$r,f_plot$lambda_K/(unique(f_count_tot$abundance[f_count_tot$species==unique_sp[s1]]/volume)*unique(f_count_tot$abundance[f_count_tot$species==unique_sp[s1]])),pch=16,col=colo[s1])
	colors=colo[s1]

	for(s2 in 1:length(unique_sp)){
		if(s1!=s2){
			colors=c(colors,colo[s2])
			ss=c(ss,s2)
			f_plot=subset(f_tot,sp1==unique_sp[s1]&sp2==unique_sp[s2])
			lines(f_plot$r,f_plot$lambda_K/(unique(f_count_tot$abundance[f_count_tot$species==unique_sp[s2]]/volume)*unique(f_count_tot$abundance[f_count_tot$species==unique_sp[s1]])),lty=1,col=colo[s2])
		}
	}
}
legend("topleft",c(paste("S=",s1," x S=",ss,sep=""),"4/3 pi r^3"),col=c(colors,"black"),bty="n",pch=c(16,NA,NA,NA),lty=c(NA,1,1,2),lwd=2)
lines(unique(f_tot$r),th_poisson,lty=2,lwd=2)
lines(unique(f_tot$r),th_poisson+thomas_cdf(unique(f_tot$r))/(200/volume),lty=3,lwd=2)
}
dev.off()
