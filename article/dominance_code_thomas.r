rm(list=ls())
graphics.off()

library(pracma) #erf function

thomas_cdf=function(r,sigma){
	a=1/(sigma*sqrt(pi))*(sigma*sqrt(pi)*erf(r/(2*sigma))-r*exp(-(r^2/(4*sigma^2))))
	return(a)
}

print(thomas_cdf(0.1,0.01))

poisson_cdf=function(r){
	a=4/3*pi*r^3
	return(a)
}

asigma=0.01
intensity_parent=c(10,100,50)
intensity_children=c(10,10,10)
intensity_total=intensity_parent*intensity_children


seq_r=seq(0.01,1,length.out=100)

dom=matrix(NA,length(intensity_parent),length(seq_r))

den_1=sum(intensity_total)*poisson_cdf(seq_r)

for (i in 1:length(intensity_parent)){
	num=intensity_total[i]*(poisson_cdf(seq_r)+thomas_cdf(seq_r,asigma)/intensity_parent[i])
	den_2=intensity_total[i]*thomas_cdf(seq_r,asigma)/intensity_parent[i]
	dom[i,]=num/(den_1+den_2)
}

pdf("dominance_thomas.pdf")
plot(0,0,xlim=range(seq_r),ylim=c(0,1),xlab="r",ylab="dominance",t="n")
lines(seq_r,dom[1,],lty=1,col="black",lwd=2.5)
lines(seq_r,dom[2,],lty=2,col="blue",lwd=2.5)
lines(seq_r,dom[3,],lty=3,col="red",lwd=2.5)
legend("bottomleft",c("sp1=100","sp2=1000","sp3=500"),lty=1:3,col=c("black","blue","red"),bty="n",lwd=2.5)
dev.off()
