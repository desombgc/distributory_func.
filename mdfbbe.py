import numpy as np
import matplotlib.pyplot as plt
import scipy.constants as cons
E1= np.linspace(0.0,0.5,100)
E2= np.linspace(0.0,2.0,100)
E3= np.linspace(0.0001,0.2,100)
Ef=1.0
kb=cons.k/(1.6e-19)
mu= lambda t: Ef*(1-(np.pi**2/12)*((kb*t)/Ef)**2)

def mb(T,E):
    return  1/np.exp(E/(kb*T))
plt.title("Maxwell-Boltzman distribution function with energy", size=12)
plt.xlabel("Energy(eV)",size= 14)
plt.ylabel("Distribution Function[f(E)]",size=14)
plt.plot(E1,mb(200,E1),marker='*',mfc='k',ms=3,ls='--',color='r',label="For T= 200K")
plt.plot(E1,mb(400,E1),marker='o',mfc='k',ms=3,ls='--',color='b',label="For T= 400K")
plt.plot(E1,mb(600,E1),marker='^',mfc='k',ms=3,ls='--',color='g',label="For T= 600K")    
plt.legend(loc="best",prop={'size':13})
plt.grid()
plt.show()

def fd(T,E):
    return 1/(np.exp((E-mu(T))/(kb*T))+1)
plt.title("Fermi-Dirac distribution function with energy", size=12)
plt.xlabel("Energy(eV)",size= 14)
plt.ylabel("Distribution Function[f(E)]",size=14)
plt.plot(E2,fd(0,E2),marker='.',mfc='k',ms=3,ls='--',color='m',label="For T= 0K")
plt.plot(E2,fd(200,E2),marker='*',mfc='k',ms=3,ls='--',color='r',label="For T= 200K")
plt.plot(E2,fd(400,E2),marker='o',mfc='k',ms=3,ls='--',color='b',label="For T= 400K")
plt.plot(E2,fd(600,E2),marker='^',mfc='k',ms=3,ls='--',color='g',label="For T= 600K")    
plt.legend(loc="best",prop={'size':13})
plt.grid()
plt.show()

def be(T,E):
    return 1/(np.exp(E/(kb*T))-1)
plt.title("Bose-Einstein distribution function with energy", size=12)
plt.ylim(0,25)
plt.xlabel("Energy(eV)",size= 14)
plt.ylabel("Distribution Function[f(E)]",size=14)
plt.plot(E3,be(200,E3),marker='*',mfc='k',ms=3,ls='--',color='r',label="For T= 200K")
plt.plot(E3,be(400,E3),marker='o',mfc='k',ms=3,ls='--',color='b',label="For T= 400K")
plt.plot(E3,be(600,E3),marker='^',mfc='k',ms=3,ls='--',color='g',label="For T= 600K")    
plt.legend(loc="best",prop={'size':13})
plt.grid()
plt.show()

