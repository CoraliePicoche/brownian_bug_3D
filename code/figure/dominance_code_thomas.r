rm(list=ls())
graphics.off()

library(scatterplot3d)

set.seed(1)

source("theoretical_functions.r")
source("utilitary_functions.r")
colo=c("red","blue","grey")
alty=c(1,2,4)

asigma=0.01
intensity_parent=c(200,100,100)
intensity_children=c(50,10,50)
intensity_total=intensity_parent*intensity_children


pdf("example_Thomas_distribution.pdf",width=5,height=5)

layout(matrix(c(1,1,2,2,3,3,4,4), c(2,4),byrow=T))

lim=1
spl=scatterplot3d(0,0,0,type="n",xlim=c(0,lim),ylim=c(0,lim),zlim=c(0,lim),xlab="x",ylab="",zlab="z",cex.axis=0.6,mar=c(2.5,2.5,2,1),cex.lab=0.7)
text(x = 8.5, y = 0.5, "y", srt = 45,cex.lab=0.7)
mtext("a",side=3,line=0.1,font=2,at=0.1,cex=0.7)

for(i in 1:length(intensity_parent)){
	x=c()
	y=c()
	z=c()
	m_p=rpois(1,intensity_parent[i])
	for (j in 1:m_p){
		xtmp=runif(1,0,1)
		ytmp=runif(1,0,1)
		ztmp=runif(1,0,1)
		x=c(x,xtmp)
		y=c(y,ytmp)
		z=c(z,ztmp)
		m_c=rpois(1,intensity_children[i])
		for (c in 1:m_c){
			x_bis=rnorm(1,xtmp,asigma)
			y_bis=rnorm(1,ytmp,asigma)
			z_bis=rnorm(1,ztmp,asigma)
			x=c(x,x_bis)
			y=c(y,y_bis)
			z=c(z,z_bis)
		}
	}
	spl$points3d(x, y, z, col=colo[i],pch=16,cex=0.5)
}

pow_r=seq(-4,0,length.out=100)
seq_r=10^pow_r

par(mar=c(4,4,4,4))
###########PCF
plot(0,0,xlim=range(seq_r),ylim=c(0,250),xlab="r",ylab="g",t="n",log="x",axes=F)
axis(1, at=log10Tck('x','major'), tcl= 0.5) # bottom
axis(1, at=log10Tck('x','minor'), tcl= 0.1, labels=NA) # bottom
axis(2)
box()
mtext("b",side=3,line=2.1,font=2,at=8e-5,cex=0.7)


for(i in 1:length(intensity_parent)){
	th_thomas=1+exp(-seq_r^2/(4*asigma^2))*(4*pi*asigma^2)^(-3/2)*1/intensity_parent[i]
	lines(seq_r,th_thomas,col=colo[i],lty=alty[i],lwd=2)
}

###########Ripley's functions
plot(0.1,0.1,xlim=range(seq_r),ylim=c(1e-7,5),xlab="r",ylab="K",t="n",log="xy",axes=F)
axis(1, at=log10Tck('x','major'), tcl= 0.5) # bottom
axis(1, at=log10Tck('x','minor'), tcl= 0.1, labels=NA) # bottom
axis(2, at=log10Tck('y','major'), tcl= 0.5)
axis(2, at=log10Tck('y','minor'), tcl= 0.1, labels=NA)
box()
mtext("c",side=3,line=2.1,font=2,at=8e-5,cex=0.7)


for(i in 1:length(intensity_parent)){
	th_thomas=thomas_cdf(seq_r,asigma,intensity_parent[i])
	print(range(th_thomas))
	lines(seq_r,th_thomas,col=colo[i],lwd=2,lty=alty[i])
}

###########Dominance
plot(0,0,xlim=range(seq_r),ylim=c(0,1),xlab="r",ylab="D",t="n",log="x",axes=F)
axis(1, at=log10Tck('x','major'), tcl= 0.5) # bottom
axis(1, at=log10Tck('x','minor'), tcl= 0.1, labels=NA) # bottom
axis(2)
box()
mtext("d",side=3,line=2.1,font=2,at=8e-5,cex=0.7)

th_poisson=4/3*pi*seq_r^3
den_1=sum(intensity_total)*th_poisson

for(i in 1:length(intensity_parent)){
	th_thomas=thomas_cdf(seq_r,asigma,intensity_parent[i])
	num=intensity_total[i]*th_thomas
	th_dom=num/(den_1+num-intensity_total[i]*th_poisson)
	lines(seq_r,th_dom,lty=alty[i],col=colo[i],lwd=2)
}

legend("bottomleft",c(paste("Sp1=",round(100*intensity_total[1]/sum(intensity_total)),"%"),paste("Sp2=",round(100*intensity_total[2]/sum(intensity_total)),"%"),paste("Sp3=",round(100*intensity_total[3]/sum(intensity_total)),"%")),col=colo,lty=alty,lwd=2,bty="n")

dev.off()
