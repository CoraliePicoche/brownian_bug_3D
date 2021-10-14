rm(list=ls())
graphics.off()

vit=0
V=0.5
x=0.5
y=0.5
z=0.5

####3D
for(i in 1:1000000){
	x_init=x
	y_init=y
	z_init=z
	x=x_init+V*cos(2*pi*y_init+runif(1,0,2*pi))
	y=y_init+V*cos(2*pi*z_init+runif(1,0,2*pi))
	z=z_init+V*cos(2*pi*x+runif(1,0,2*pi))
	vit_tmp=(x-x_init)^2+(y-y_init)^2+(z-z_init)^2
#	vit=c(vit,sqrt(vit_tmp))
	vit=vit+sqrt(vit_tmp)
}
vit=vit/10^6 #vit=V*1.194    V=vit*0.837 
print(vit/V)


vit=0
#Let's try something in 2D
for(i in 1:1000000){
        x_init=x
        y_init=y
        x=x_init+0.5*V*cos(2*pi*y_init+runif(1,0,2*pi))
        y=y_init+0.5*V*cos(2*pi*x+runif(1,0,2*pi))
        vit_tmp=(x-x_init)^2+(y-y_init)^2
#       vit=c(vit,sqrt(vit_tmp))
        vit=vit+sqrt(vit_tmp)
}
vit=vit/10^6 ##vit=V*0.958    V=vit*1.04 (V=vit*pi/3 ??????)
print(vit/V)
