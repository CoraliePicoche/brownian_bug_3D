rm(list=ls())
graphics.off()

library(scatterplot3d)

source("theoretical_functions.r")
colo=c("red","blue","grey")
alty=c(1,2,4)

asigma=0.01
intensity_parent=c(200,100,100)
intensity_children=c(50,10,50)
intensity_total=intensity_parent*intensity_children


pdf("example_Thomas_distribution.pdf",width=5,height=5)
par(cex=1.0,oma=c(2,2,2,2),mar=c(2,2,2,2),xpd=T)

#par(mfrow=c(2,2))
layout(matrix(c(1,1,2,2,3,3,4,4), c(2,4),byrow=T))

lim=1
spl=scatterplot3d(0,0,0,type="n",xlim=c(0,lim),ylim=c(0,lim),zlim=c(0,lim),xlab="x",ylab="",zlab="z",main="Thomas process",cex.axis=0.6,mar=c(2.5,2.5,2,1),cex.lab=0.7)
text(x = 8.5, y = 0.5, "y", srt = 45,cex.lab=0.7)

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

pow_r=seq(-3,0,length.out=100)
seq_r=10^pow_r

###########PCF
plot(0,0,xlim=range(seq_r),ylim=c(0,250),xlab="r",ylab="g(r)",t="n",log="x")

for(i in 1:length(intensity_parent)){
	th_thomas=1+exp(-seq_r^2/(4*asigma^2))*(4*pi*asigma^2)^(-3/2)*1/intensity_parent[i]
	lines(seq_r,th_thomas,col=colo[i],lty=alty[i],lwd=2)
}

###########Ripley's functions
plot(0.1,0.1,xlim=range(seq_r),ylim=c(1e-7,5),xlab="r",ylab="K(r)",t="n",log="xy")

for(i in 1:length(intensity_parent)){
	th_thomas=thomas_cdf(seq_r,asigma,intensity_parent[i])
	print(range(th_thomas))
	lines(seq_r,th_thomas,col=colo[i],lwd=2,lty=alty[i])
}

###########Dominance
plot(0,0,xlim=range(seq_r),ylim=c(0,1),xlab="r",ylab="dominance",t="n",log="x")

th_poisson=4/3*pi*seq_r^3
den_1=sum(intensity_total)*th_poisson

for(i in 1:length(intensity_parent)){
	th_thomas=thomas_cdf(seq_r,asigma,intensity_parent[i])
	num=intensity_total[i]*th_thomas
	th_dom=num/(den_1+num-intensity_total[i]*th_poisson)
	lines(seq_r,th_dom,lty=alty[i],col=colo[i],lwd=2)
}
dev.off()
