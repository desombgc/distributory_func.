import numpy as np
import matplotlib.pyplot as plt
import scipy.constants as cons
from scipy.integrate import quad

R=cons.R
n=100
T=np.linspace(1,700,n)
sh=3*R

def cv(T):
    return sh*(T/T)

def cve(Te,T):
    return 3*R*(Te/T)**2*(np.exp(Te/T)/(np.exp(Te/T)-1)**2)

def cvd(Td,T):
 shd=[]
 p=lambda x: ((x**4)*np.exp(x))/((np.exp(x)-1)**2)
 for i in T:  
       res=quad(p,0,Td/i)[0]
       shd.append(9*R*(i/Td)**3*res)
 return shd
   
plt.title("Dulong-Pettit Law", size= 16)
plt.xlabel("Temperature($T$) in K", size=14)
plt.ylabel("$C$$_V$ (Joule/mol.K)", size=14)
plt.plot(T,cv(T),marker='*',mfc='k',ms=3,ls='--',color='purple',label="Dulong-Pettit Law")
plt.legend(loc="best",prop={'size':13})
plt.grid()
plt.show()
plt.title("Einstein's law of specific heat at constant volume", size=16)
plt.xlabel("Temperature($T$) in K", size=14)
plt.ylabel("$C$$_V$ (Joule/mol.K)", size=14)
plt.plot(T,cve(160,T),marker='*',mfc='k',ms=3,ls='--',color='r',label="$\\theta$$_E$=160K(Silver)")
plt.plot(T,cve(240,T),marker='^',mfc='k',ms=3,ls='--',color='g',label="$\\theta$$_E$=240K(Copper)")
plt.plot(T,cve(290,T),marker='o',mfc='k',ms=3,ls='--',color='b',label="$\\theta$$_E$=290K(Aluminium)")
plt.legend(loc="best",prop={'size':13})
plt.grid()
plt.show()
plt.title("Debye law of specific heat at constant volume", size=16)
plt.xlabel("Temperature($T$) in K", size=14)
plt.ylabel("$C$$_V$ (Joule/mol.K)", size=14)
plt.plot(T,cvd(215,T),marker='*',mfc='k',ms=3,ls='--',color='r',label="$\\theta$$_D$=215K(Silver)")
plt.plot(T,cvd(348,T),marker='^',mfc='k',ms=3,ls='--',color='b',label="$\\theta$$_D$=348K(Copper)")
plt.plot(T,cvd(428,T),marker='o',mfc='k',ms=3,ls='--',color='g',label="$\\theta$$_D$=428K(Aluminium)")
plt.legend(loc="best",prop={'size':13})
plt.grid()
plt.show()
