import numpy as np
import matplotlib.pyplot as plt
import scipy as sc
import matplotlib.cm as cm

def BBM_cdf(r,gamma,D,Conc,lambda_val,t=np.nan):
        res=np.nan
        A=4/3*np.pi*r**3
        if(gamma>0):
                B=(r**2)/(6*D)
                C=(np.sqrt(gamma)/(6*np.sqrt(2)*(D**1.5)))*(r**3)*np.arctan(r*np.sqrt(gamma/(2*D)))
                Dd=np.log((gamma*r**2)/(2*D)+1)/(6*gamma)
                E=(np.sqrt(gamma)*np.pi*r**3)/(12*np.sqrt(2)*D**1.5)
                res=A+(2*lambda_val/Conc)*(B+C+Dd-E)
        else:
                B=(r**2)/2
                C=1/2*sc.erf(r/(np.sqrt(8*D*t)))*(r**2-4*D*t)
                Dd=np.sqrt(2*D*t)*r/np.sqrt(np.pi)*np.exp(-r**2/(8*D*t))
                res=A+(lambda_val/(Conc*D))*(B-C-Dd)
        
        return res

#Parameters Diatoms
a_lambda=1*0.0002
a_C=10000/1000
a_r=25*10**(-6)

#Common parameters
a_D=8.314*293/(6.0225*10**23)*1/(6.0*np.pi*a_r*10**(-3))
n_species=3
seq_r=np.linspace(10**(-4),10**(-3),100)
seq_r=np.hstack([seq_r,np.linspace(10**(-3),10**(-2),100)])
seq_r=np.hstack([seq_r,np.linspace(10**(-2),10**(-1),100)])
seq_r=np.hstack([seq_r,np.linspace(10**(-1),10**(0),100)])

a_gamma=0.5
val_K=np.empty((1,400))

for i in range(len(seq_r)):
    val_K[0,i]=BBM_cdf(seq_r[i],a_gamma,a_D,a_C,a_lambda)

#print(len(seq_r))
print(seq_r[0:5])
print(seq_r[394:400])
#print(val_K)
fig=plt.figure()
plt.plot(seq_r,val_K[0,:])
plt.xlabel("r",fontsize="x-small")
plt.ylabel("K",fontsize="x-small")
plt.savefig('K_func_BBM_python.pdf',bbox_inches='tight')

