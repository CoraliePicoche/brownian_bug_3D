rm(list=ls())
graphics.off()

set.seed(42)

velocity=0
U=0.5
Lmax=10
max_v=0

####3D
#Initial position of the particle
for (s in 1:1000){
	x=runif(1,0,Lmax)
	y=runif(1,0,Lmax)
	z=runif(1,0,Lmax)
	for(i in 1:1000){
		x_init=x
		y_init=y
		z_init=z
		x=x_init+(1/2)*U*cos(2*pi*y_init+runif(1,0,2*pi))
		y=y_init+(1/2)*U*cos(2*pi*z_init+runif(1,0,2*pi))
		z=z_init+(1/2)*U*cos(2*pi*x+runif(1,0,2*pi))
		vit_tmp=(x-x_init)^2+(y-y_init)^2+(z-z_init)^2
		velocity=velocity+sqrt(vit_tmp)
	}
}
velocity=velocity/10^6 
print(velocity/U)

####2D
velocity=0
for (s in 1:1000){
	x=runif(1,0,Lmax)
	y=runif(1,0,Lmax)
	z=runif(1,0,Lmax)
	for(i in 1:1000){
        	x_init=x
	        y_init=y
        	x=x_init+(1/2)*U*cos(2*pi*y_init+runif(1,0,2*pi))
	        y=y_init+(1/2)*U*cos(2*pi*x+runif(1,0,2*pi))
        	vit_tmp=(x-x_init)^2+(y-y_init)^2
	        velocity=velocity+sqrt(vit_tmp)
		if(vit_tmp>max_v) max_v=vit_tmp
	}
}
velocity=velocity/10^6
print(velocity/U)
