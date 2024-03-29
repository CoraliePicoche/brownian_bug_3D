rm(list=ls())
graphics.off()

library(scatterplot3d)

#To have different shades for a given species
#paletteFunc <- c(colorRampPalette(c('lightblue', 'darkblue')),colorRampPalette(c('orange', 'darkred')),colorRampPalette(c('lightgreen', 'darkgreen')),colorRampPalette(c('plum2', 'darkorchid4')))


png("spatial_distribution_micro.png",width=1000)
par(mfrow=c(1,2),cex=1.5,xpd=T)

f=read.table("../simulation/Spatial_distribution_0.txt",header=F,sep=";",dec='.')
colnames(f)=c("id","x","y","z","yfirst","first_parent","species")
unique_sp=unique(f$species)
nb_sp=length(unique_sp)

palette=c("red","blue","grey")

spl=scatterplot3d(0,0,0,type="n",xlim=c(0,10),ylim=c(0,10),zlim=c(0,10),xlab="x",ylab="",zlab="z",main="Advection",cex.lab=1.5,cex.axis=1.25)
text(x = 8, y = 0.75, "y", srt = 45,cex.lab=0.7)
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
spl=scatterplot3d(0,0,0,type="n",xlim=c(0,10),ylim=c(0,10),zlim=c(0,10),xlab="x",ylab="",zlab="z",main="No advection",cex.lab=1.5,cex.axis=1.25)
text(x = 8, y = 0.75, "y", srt = 45,cex.lab=0.7)
for (s in 1:nb_sp){
        ftmp=subset(f,species==unique_sp[s])
	
	id=sample(1:nrow(ftmp),1000)

	x=ftmp$x[id]
	y=ftmp$y[id]
	z=ftmp$z[id]
	
        spl$points3d(x, y, z, col=palette[s],pch=16,cex=0.5)
}

dev.off()

png("spatial_distribution_nano.png",width=1000)
par(mfrow=c(1,2),cex=1.5,xpd=T)

f=read.table("../simulation/Spatial_distribution_2.txt",header=F,sep=";",dec='.')
colnames(f)=c("id","x","y","z","yfirst","first_parent","species")
unique_sp=unique(f$species)
nb_sp=length(unique_sp)

spl=scatterplot3d(0,0,0,type="n",xlim=c(0,2),ylim=c(0,2),zlim=c(0,2),xlab="x",ylab="",zlab="z",main="Advection",cex.lab=1.3,cex.axis=1.1)
text(x = 6.5, y = .6, "y", srt = 45,cex=0.85)
for (s in 1:nb_sp){
        ftmp=subset(f,species==unique_sp[s] & x<=2 & y<=2 & z<=2)

        id=sample(1:nrow(ftmp),1000)

        x=ftmp$x[id]
        y=ftmp$y[id]
        z=ftmp$z[id]

        spl$points3d(x, y, z, col=palette[s],pch=16,cex=0.5)
}

f=read.table("../simulation/Spatial_distribution_3.txt",header=F,sep=";",dec='.')
colnames(f)=c("id","x","y","z","yfirst","first_parent","species")
unique_sp=unique(f$species)
nb_sp=length(unique_sp)
spl=scatterplot3d(0,0,0,type="n",xlim=c(0,2),ylim=c(0,2),zlim=c(0,2),xlab="x",ylab="",zlab="z",main="No advection",cex.lab=1.3,cex.axis=1.1)
text(x = 6.5, y = .6, "y", srt = 45,cex=0.85)
for (s in 1:nb_sp){
        ftmp=subset(f,species==unique_sp[s]& x<=2 & y<=2 & z<=2)

        id=sample(1:nrow(ftmp),1000)

        x=ftmp$x[id]
        y=ftmp$y[id]
        z=ftmp$z[id]

        spl$points3d(x, y, z, col=palette[s],pch=16,cex=0.5)
}

dev.off()

###Now with a zoom on micro

lim_min=0
lim_max=lim_min+0.5

while(lim_min<lim_max){

pdf(paste("spatial_distribution_zoom_micro",lim_min,".pdf",sep=""),width=5,height=5)
par(cex=1.0,oma=c(0,0,0,0),mar=c(0,0,0,0),xpd=T)

layout(matrix(c(1,1,2,2,3,3,4,4), c(2,4),byrow=T))

f=read.table("../simulation/Spatial_distribution_0.txt",header=F,sep=";",dec='.')
colnames(f)=c("id","x","y","z","yfirst","first_parent","species")
unique_sp=unique(f$species)
nb_sp=length(unique_sp)

spl=scatterplot3d(lim_min,lim_min,lim_min,type="n",xlim=c(lim_min,lim_min+2),ylim=c(lim_min,lim_min+2),zlim=c(lim_min,lim_min+2),xlab="x",ylab="",zlab="z",main="Advection",cex.axis=0.6,mar=c(2.5,2.5,2,1),cex.lab=0.7)
text(x = 6.5, y = 0.6, "y", srt = 45,cex.lab=0.7)
for (s in 1:nb_sp){
        ftmp=subset(f,species==unique_sp[s] & x>=lim_min & x<=(lim_min+2) & y>=lim_min & y<=(lim_min+2) & z>=lim_min & z<=(lim_min+2))
        x=ftmp$x
        y=ftmp$y
        z=ftmp$z

        spl$points3d(x, y, z, col=palette[s],pch=16,cex=0.5)
}

par(mar=c(4,4,4,4))
plot(0,0,type="n",xlim=c(lim_min,lim_min+2),ylim=c(lim_min,lim_min+2),xlab="x",ylab="y")
for (s in 1:nb_sp){
        ftmp=subset(f,species==unique_sp[s] & x>=lim_min & x<=(lim_min+2) & y>=lim_min & y<=(lim_min+2) & z>=lim_min & z<=(lim_min+2))

        x=ftmp$x
        y=ftmp$y

        points(x, y, col=palette[s],pch=16,cex=0.5)
}


f=read.table("../simulation/Spatial_distribution_1.txt",header=F,sep=";",dec='.')
colnames(f)=c("id","x","y","z","yfirst","first_parent","species")
unique_sp=unique(f$species)
nb_sp=length(unique_sp)

spl=scatterplot3d(lim_min,lim_min,lim_min,type="n",xlim=c(lim_min,lim_min+2),ylim=c(lim_min,lim_min+2),zlim=c(lim_min,lim_min+2),xlab="x",ylab="",zlab="z",main="No advection",cex.axis=0.6,mar=c(2.5,2.5,2,1),cex.lab=0.7)
text(x = 6.5, y = 0.6, "y", srt = 45,cex.lab=0.7)
for (s in 1:nb_sp){
        ftmp=subset(f,species==unique_sp[s] & x>=lim_min & x<=(lim_min+2) & y>=lim_min & y<=(lim_min+2) & z>=lim_min & z<=(lim_min+2))

        x=ftmp$x
        y=ftmp$y
        z=ftmp$z

        spl$points3d(x, y, z, col=palette[s],pch=16,cex=0.5)
}

par(mar=c(4,4,4,4))
plot(0,0,type="n",xlim=c(lim_min,lim_min+2),ylim=c(lim_min,lim_min+2),xlab="x",ylab="y")
for (s in 1:nb_sp){
        ftmp=subset(f,species==unique_sp[s] & x>=lim_min & x<=(lim_min+2) & y>=lim_min & y<=(lim_min+2) & z>=lim_min & z<=(lim_min+2))

        x=ftmp$x
        y=ftmp$y

        points(x, y, col=palette[s],pch=16,cex=0.5)

}

dev.off()
lim_min=lim_min+1
}

############Zoom on nano

lim=0.5

pdf("spatial_distribution_zoom_nano.pdf",width=5,height=5)
par(cex=1.0,oma=c(0,0,0,0),mar=c(0,0,0,0),xpd=T)

#par(mfrow=c(2,2))
layout(matrix(c(1,1,2,2,3,3,4,4), c(2,4),byrow=T))

f=read.table("../simulation/Spatial_distribution_2.txt",header=F,sep=";",dec='.')
colnames(f)=c("id","x","y","z","yfirst","first_parent","species")
unique_sp=unique(f$species)
nb_sp=length(unique_sp)

spl=scatterplot3d(0,0,0,type="n",xlim=c(0,lim),ylim=c(0,lim),zlim=c(0,lim),xlab="x",ylab="",zlab="z",main="Advection",cex.axis=0.6,mar=c(2.5,2.5,2,1),cex.lab=0.7)
text(x = 8.5, y = 0.5, "y", srt = 45,cex.lab=0.7)
for (s in 1:nb_sp){
        ftmp=subset(f,species==unique_sp[s] & x<=lim & y<=lim & z<=lim)

        x=ftmp$x
        y=ftmp$y
        z=ftmp$z

        spl$points3d(x, y, z, col=palette[s],pch=16,cex=0.5)
}

par(mar=c(4,4,4,4))
plot(0,0,type="n",xlim=c(0,lim),ylim=c(0,lim),xlab="x",ylab="y",cex.lab=1,cex.axis=1)
for (s in 1:nb_sp){
        ftmp=subset(f,species==unique_sp[s] & x<=lim & y<=lim & z<=lim)

        x=ftmp$x
        y=ftmp$y

        points(x, y, col=palette[s],pch=16,cex=0.5)
}

f=read.table("../simulation/Spatial_distribution_3.txt",header=F,sep=";",dec='.')
colnames(f)=c("id","x","y","z","yfirst","first_parent","species")
unique_sp=unique(f$species)
nb_sp=length(unique_sp)

spl=scatterplot3d(0,0,0,type="n",xlim=c(0,lim),ylim=c(0,lim),zlim=c(0,lim),xlab="x",ylab="",zlab="z",main="No advection",cex.axis=0.6,mar=c(2.3,2.5,2.2,1),cex.lab=0.7)
text(x = 8.5, y = 0.5, "y", srt = 45,cex.lab=0.7)
for (s in 1:nb_sp){
        ftmp=subset(f,species==unique_sp[s] & x<=lim & y<=lim & z<=lim)

        x=ftmp$x
        y=ftmp$y
        z=ftmp$z

        spl$points3d(x, y, z, col=palette[s],pch=16,cex=0.5)
}

par(mar=c(4,4,4,4))
plot(0,0,type="n",xlim=c(0,lim),ylim=c(0,lim),xlab="x",ylab="y")
for (s in 1:nb_sp){
        ftmp=subset(f,species==unique_sp[s] & x<=lim & y<=lim & z<=lim)

        x=ftmp$x
        y=ftmp$y

        points(x, y, col=palette[s],pch=16,cex=0.5)
}

dev.off()
