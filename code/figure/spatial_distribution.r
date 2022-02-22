rm(list=ls())
graphics.off()

library(scatterplot3d)

#To have different shades for a given species
#paletteFunc <- c(colorRampPalette(c('lightblue', 'darkblue')),colorRampPalette(c('orange', 'darkred')),colorRampPalette(c('lightgreen', 'darkgreen')),colorRampPalette(c('plum2', 'darkorchid4')))


png("spatial_distribution.png",width=1000)
par(mfrow=c(1,2),cex=1.5)

f=read.table("../simulation/Spatial_distribution_0.txt",header=F,sep=";",dec='.')
colnames(f)=c("id","x","y","z","yfirst","first_parent","species")
unique_sp=unique(f$species)
nb_sp=length(unique_sp)

palette=c("red","blue","grey")

spl=scatterplot3d(0,0,0,type="n",xlim=c(0,10),ylim=c(0,10),zlim=c(0,10),xlab="x",ylab="",zlab="z",main="Advection")
for (s in 1:nb_sp){
	ftmp=subset(f,species==unique_sp[s])
	
	id=sample(1:nrow(ftmp),1000)

	x=ftmp$x[id]
	y=ftmp$y[id]
	z=ftmp$z[id]

	spl$points3d(x, y, z, col=palette[s],pch=16,cex=0.5)
}

f=read.table("../simulation/Spatial_distribution_1.txt",header=F,sep=";",dec='.')
colnames(f)=c("id","x","y","z","yfirst","first_parent","species")
unique_sp=unique(f$species)
nb_sp=length(unique_sp)
spl=scatterplot3d(0,0,0,type="n",xlim=c(0,10),ylim=c(0,10),zlim=c(0,10),xlab="x",ylab="y",zlab="z",main="No advection")
for (s in 1:nb_sp){
        ftmp=subset(f,species==unique_sp[s])
	
	id=sample(1:nrow(ftmp),1000)

	x=ftmp$x[id]
	y=ftmp$y[id]
	z=ftmp$z[id]
	
        spl$points3d(x, y, z, col=palette[s],pch=16,cex=0.5)
}


dev.off()
