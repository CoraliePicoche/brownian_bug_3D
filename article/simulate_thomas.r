rm(list=ls())
graphics.off()

set.seed(22)

library(pracma)

thomas_cdf=function(r,sigma,d){
	if(d==2){
		a=(1-exp(r^2/(4*sigma)^2))
	}else if(d==3){
        	a=2/(sigma*sqrt(pi))*(erf(r/(2*sigma))-r/(sigma*sqrt(pi))*exp(-(r^2/(4*sigma^2))))
	}
	return(a)
}

poisson_cdf=function(r,d){
	if(d==2){
		a=pi*r^2
	}else if(d==3){
        	a=4/3*pi*r^3
	}
        return(a)
}

asigma=0.01
D=2

particle_matrix=matrix(NA,,4)
colnames(particle_matrix)=c("species","x","y","z")

intensity_parent=c(10,100,50)
sp=length(intensity_parent)
intensity_children=c(10,10,10)
actual_intensity_children=rep(0,sp)

dist_mat=c()
dist_mat_inter=c()

for(s in 1:sp){
	for(p in 1:intensity_parent[s]){
		x_p=runif(1,0,1)
		y_p=runif(1,0,1)
		z_p=runif(1,0,1)
		nb_child=rpois(1,intensity_children[s])
		actual_intensity_children[s]=actual_intensity_children[s]+nb_child
		vec_p=c(s,x_p,y_p,z_p)
		particle_matrix=rbind(particle_matrix,vec_p)
		for (c in 1:nb_child){
			x_c=rnorm(1,x_p,asigma)
			y_c=rnorm(1,y_p,asigma)
			z_c=rnorm(1,z_p,asigma)
			vec_c=c(s,x_c,y_c,z_c)
			particle_matrix=rbind(particle_matrix,vec_c)
		}
	}
}
particle_matrix=particle_matrix[2:nrow(particle_matrix),]

seq_r=seq(0.01,1,length.out=10)
K_matrix=array(0,dim=list(length(seq_r),sp,sp))

par(mfrow=c(1,2))
colo=c("black","blue","red")
plot(0,0,xlim=c(0,1),ylim=c(0,1),xlab="",ylab="")
for(s1 in 1:sp){
	s2=s1
	while(s2<=sp){
		tt_1=subset(particle_matrix,particle_matrix[,"species"]==s1)
		tt_2=subset(particle_matrix,particle_matrix[,"species"]==s2)
		if(s2==s1){
			points(tt_1[,"x"],tt_1[,"y"],col=colo[s1])
		}
		for(i in 1:nrow(tt_1)){
			for(j in 1:nrow(tt_2)){
				if(j!=i){
					dist=(tt_1[i,"x"]-tt_2[j,"x"])^2+(tt_1[i,"y"]-tt_2[j,"y"])^2
					if(D==3){
						dist=dist+(tt_1[i,"z"]-tt_2[j,"z"])^2
					}
					r=1
					while (r<length(seq_r) & dist>seq_r[r]^2){
						r=r+1
					}
					K_matrix[r:length(seq_r),s1,s2]=K_matrix[r:length(seq_r),s1,s2]+1
					#K_matrix[r,s1,s2]=K_matrix[r,s1,s2]+1
				}			
			}	
		}
	s2=s2+1
	}
}

plot(0,0,xlim=range(seq_r),ylim=c(0,4/3*pi),xlab="r",ylab="K",t="n")
lines(seq_r,poisson_cdf(seq_r,D))
#lines(seq_r,intensity_parent[1]*intensity_children[1]*(poisson_cdf(seq_r)+thomas_cdf(seq_r,asigma)/intensity_parent[1]),col="green")
lines(seq_r,(poisson_cdf(seq_r,D)+thomas_cdf(seq_r,asigma,D)/intensity_parent[1]),col="green")
points(seq_r,cumsum(K_matrix[,1,1])/((intensity_parent[1]*intensity_children[1])^2),col="black",pch=16)
points(seq_r,cumsum(K_matrix[,1,2])/((intensity_parent[1]*intensity_children[1])*intensity_parent[2]*intensity_children[2]),col="red",pch=16)
stop()

plot(0,0,xlim=range(seq_r),ylim=c(0,4/3*pi),xlab="r",ylab="K",t="n")
lines(seq_r,4/3*pi*seq_r^3)
points(seq_r,cumsum(K_matrix[,2,2])/((intensity_parent[2]*intensity_children[2])^2),col="black",pch=16)
points(seq_r,cumsum(K_matrix[,2,3])/((intensity_parent[3]*intensity_children[3])*intensity_parent[2]*intensity_children[2]),col="red",pch=16)

plot(0,0,xlim=range(seq_r),ylim=c(0,4/3*pi),xlab="r",ylab="K",t="n")
lines(seq_r,4/3*pi*seq_r^3)
points(seq_r,cumsum(K_matrix[,3,3])/((intensity_parent[3]*intensity_children[3])^2),col="black",pch=16)
points(seq_r,cumsum(K_matrix[,1,3])/((intensity_parent[1]*intensity_children[1])*intensity_parent[3]*intensity_children[3]),col="red",pch=16)
