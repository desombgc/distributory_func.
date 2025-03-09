import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import simps
import scipy.constants as cons
from scipy.optimize import fmin

pi=cons.pi
c=cons.c
kb=cons.k
h=cons.h
T0=4000;T1=8000;T2=300
nu0=np.linspace(1,15*kb*T0/h,200)
nu1=np.linspace(1,15*kb*T1/h,200)
nu2=np.linspace(1,15*kb*T2/h,200)

def urj(nu,T):
    return  (8*pi*nu**2*kb*T)/(c**3)
def upk(nu,T):
    return  (8*pi*nu**3*h)/((c**3)*(np.exp((h*nu)/(kb*T))-1))

##PLanck's law and rayleigh-jeans law for blackbody radiation for 4000K and 8000K temp
plt.ylim(0,6e-16)
plt.title("Planck's law and Rayliegh-Jeans law at T=4000K", size=15)
plt.grid()
plt.plot(nu0,urj(nu0,T0),marker='^',mfc='r',ms='3',ls='-.',label="Rayleigh-Jeans Law")
plt.plot(nu0,upk(nu0,T0),marker='*',mfc='k',ms='3',ls='--',label="Planck's Law")
plt.xlabel("Frequency($\\nu$)",size=15)
plt.ylabel("Energy density(u$_\\nu$)",size=15)
plt.legend(loc="best",prop={'size':10})
plt.show()

plt.ylim(0,5e-15)
plt.title("Planck's law and Rayliegh-Jeans law at T=8000K", size=15)
plt.grid()
plt.plot(nu1,urj(nu1,T1),marker='^',mfc='r',ms='3',ls='-.',label="Rayleigh-Jeans Law")
plt.plot(nu1,upk(nu1,T1),marker='*',mfc='k',ms='3',ls='--',label="Planck's Law")
plt.xlabel("Frequency($\\nu$)",size=15)
plt.ylabel("Energy density(u$_\\nu$)",size=15)
plt.legend(loc="best",prop={'size':10})
plt.show()
##frequency at which the energy density is maximum
x0=(h*nu0)/(kb*T0)
x1=(h*nu1)/(kb*T1)
def p(x0):
    return ((x0)**3)*(np.exp(-x0))/(1-np.exp(-x0))
def q(x1):
    return ((x1)**3)*(np.exp(-x1))/(1-np.exp(-x1))
max1=fmin(lambda x0: -p(x0),2)[0]
max2=fmin(lambda x1: -p(x1),2)[0]

M1=kb*T0*max1/h
M2=kb*T1*max2/h

print("Frequency at which the energy density is maximum at T=4000 is = ",M1/1e12,'Terahertz.')
print("Frequency at which the energy density is maximum at T=8000 is = ",M2/1e12,'Terahertz.')
###C.MINIMUM frequency at where the difference between energy densities are more than 10% at T=300K
dU=(urj(nu2,T2)-upk(nu2,T2))/upk(nu2,T2)
for i in (dU):
    if i >= 0.1:
        a=list(dU).index(i)
        print("Minimum frequency where the difference between the energy densities are >= 10% is at T=300 = ",nu2[a]/1e12,'Terahertz.')
        break
    
###PLotting the visible energy as a function of temperature for Rayleigh-Jeans Law & Plancks Law
T=np.linspace(100,400,1000)
Erj,Ep=[],[]
for i in (T):
     nu=np.linspace(4e14, 7.5e14,1000)  
     Urj=(8*pi*nu**2*kb*i)/(c**3)
     Upk=(8*pi*nu**3*h)/((c**3)*(np.exp((h*nu)/(kb*i))-1)) 
     E1=simps(Urj,nu);E2=simps(Upk,nu)
     Erj.append(E1);Ep.append(E2)     

plt.title("Total energy in the visible spectrum (4e14-7.5e14 Hz),T=100-400K",size=14)
plt.xlabel("Temperature",size=15)
plt.ylabel("Energy(T)",size=15)
plt.plot(T,Erj,marker='*',mfc='g',ms='3',color='red',ls='-.',label='Rayleigh-Jeans Law')
plt.legend(loc="best",prop={'size':13})
plt.grid()
plt.show()
    
plt.title("Total energy in the visible spectrum (4e14-7.5e14 Hz),T=100-400K",size=14)
plt.xlabel("Temperature",size=15)
plt.ylabel("Energy(T)",size=15)
plt.plot(T,Ep,marker='^',mfc='b',ms='3',ls='-.',color='blue',label='Plancks Law')
plt.legend(loc="best",prop={'size':13})
plt.grid()
plt.show()  
              