rm(list=ls())
graphics.off()

asigma=0.01

particle_matrix=matrix(NA,,4)
colnames(particle_matrix)=c("species","x","y","z")

intensity_parent=c(10,100,50)
intensity_children=c(10,10,10)

sp=length(intensity_parent)

for(s in 1:sp){
	for(p in 1:intensity_parent[s]){
		x_p=runif(1,0,1)
		y_p=runif(1,0,1)
		z_p=runif(1,0,1)
		nb_child=rpois(1,intensity_children[s])
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
		tt=subset(particle_matrix,particle_matrix[,"species"]==s1 | particle_matrix[,"species"]==s2)
		if(s2==s1){
			points(tt[,"x"],tt[,"y"],col=colo[s1])
		}
		for(i in 1:nrow(tt)){
			for(j in 1:nrow(tt)){
				if(j!=i){
					dist=(tt[i,"x"]-tt[j,"x"])^2+(tt[i,"y"]-tt[j,"y"])^2+(tt[i,"z"]-tt[j,"z"])^2
					for (r in 1:length(seq_r)){
						if(dist<seq_r[r]^2){
							K_matrix[r,s1,s2]=K_matrix[r,s1,s2]+1
						}
					}
				}			
			}	
		}
	s2=s2+1
	}
}

plot(0,0,xlim=range(seq_r),ylim=range(c(K_matrix)),xlab="r",ylab="K")
lines(seq_r,4/3*pi*seq_r^3)
points(seq_r,K_matrix[,1,1],col="black")
points(seq_r,K_matrix[,1,2]/(intensity_parent[2]*intensity_children[2]),col="grey",pch=3)

