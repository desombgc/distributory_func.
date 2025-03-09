import numpy as np
import matplotlib.pyplot as plt

Ti,Tf=0.01,80.0
EN,n=50,400             #no. of energy levels, no. of points
H=(Tf-Ti)/n
kb,e=1.0,1.0

def z(i,T):
    return np.exp((-i*e)/(kb*T))

def u(i,T):
    return (i*e)*z(i,T)
def U2(i,T):
    return z(i,T)*(i*e)**2

def Z(T):
   sumz=0.0
   for i in range(1,EN):
       sumz=sumz+z(i,T)
   return sumz

def U(T):
    sumu=0.0
    for i in range(1,EN):
        sumu=sumu+u(i,T)
    return sumu/Z(T)

def U3(T):
    sumu3=0.0
    for i in range(1,EN):
        sumu3=sumu3+U2(i,T)
    return sumu3/Z(T)
def S(T):
    return np.sqrt(U3(T)-(U(T)*U(T)))
     
def cv(T):
     return S(T)**2/T**2

def F(T):
    return -kb*T*np.log(Z(T))

def en(T):
     return (U(T)-F(T))/T

def op1(T):
    return np.exp(-e/(kb*T))/Z(T)

def op2(T):
    return np.exp((-e*2)/(kb*T))/Z(T)

par=[]
temp=[]
avg_e=[]
fluc=[]
sh=[]
free=[]
entropy=[]
Pr1=[]
Pr2=[]

for j in range(1,n):
    T=Ti+j*H
    temp.append(T)
    par.append(Z(T))
    avg_e.append(U(T))
    fluc.append(S(T))
    sh.append(cv(T))
    free.append(F(T))
    entropy.append(en(T))
    Pr1.append(op1(T))
    Pr2.append(op2(T))
    
plt.title("Internal Energy(U) vs Temperature(T) Plot",size=15)
plt.xlabel("Temperature($T$)", size=14)
plt.ylabel("Internal Energy($U$)",size=14)
plt.plot(temp,avg_e,marker='^',mfc='k',ms=3,ls='--',color='g',label='For Internal Energy($U$)')
plt.legend(loc='best',prop={'size':14}) 
plt.show()
plt.title("Energy Fluctuation(ΔE) vs Temperature(T) Plot",size=15)
plt.xlabel("Temperature($T$)", size=14)
plt.ylabel("Energy Fluctuation ($ΔE$)",size=14)
plt.plot(temp,fluc,marker='*',mfc='k',ms=3,ls='--',color='purple',label='For Energy Fluctuation($ΔE$)')
plt.legend(loc='best',prop={'size':14}) 
plt.show()
plt.title("Specific Heat(Cv) vs Temperature(T) Plot",size=15)
plt.xlabel("Temperature($T$)", size=14)
plt.ylabel("Specific Heat($C_v$)",size=14)
plt.plot(temp,sh,marker='.',mfc='r',ms=3,ls='--',color='r',label='For Specific Heat($C_v$)')
plt.legend(loc='best',prop={'size':14}) 
plt.show()    
plt.title("Free Energy(F) vs Temperature(T) Plot",size=15)
plt.xlabel("Temperature($T$)", size=14)
plt.ylabel("Free Energy($F$)",size=14)
plt.plot(temp,free,marker='+',mfc='k',ms=3,ls='--',color='magenta',label='For Free Energy($F$)')
plt.legend(loc='best',prop={'size':14}) 
plt.show()    
plt.title("Entropy(S) vs Temperature(T) Plot",size=15)
plt.xlabel("Temperature($T$)", size=14)
plt.ylabel("Entropy($S$)",size=14)
plt.plot(temp,entropy,marker='v',mfc='k',ms=3,ls='--',color='b',label='For Entropy($S$)')
plt.legend(loc='best',prop={'size':14}) 
plt.show()
plt.title("Occupation number vs Temperature(T) Plot",size=15)
plt.xlabel("Temperature($T$)", size=14)
plt.ylabel("Occupation number($Pr_n$)",size=14)  
plt.plot(temp,Pr1,label='For occupation no.($Pr1$)')
plt.plot(temp,Pr2,label='For occupation no.($Pr2$)')
plt.legend(loc='best',prop={'size':14})
plt.show()