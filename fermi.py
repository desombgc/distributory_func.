import numpy as np
from itertools import combinations
import matplotlib.pyplot as plt
Ti,Tf=0.01,50
N,m=65,3          #no. of enegy level, no. of fermions 
n=100
h=(Tf-Ti)/n
kb,gamma=1.0,1.0

ek=np.ones(N)
for j in range (N):
      ek[j]=gamma*(j+1)

energy=[]
comb=combinations(ek,m)
for k in list(comb):
    energy.append(sum(np.array(k)))

E=np.array(energy)
#print(E)
def z(i,T):
    return np.exp(E[i]/-(kb*T))

def u(i,T):
    return (E[i])*z(i,T)

def u2(i,T):
    return z(i,T)*(E[i])**2

def Z(T):
    sumz=0.0
    for i in range(len(E)):
        sumz=sumz+z(i,T)
    return sumz

def U(T):
    sumu=0.0
    for i in range(len(E)):
        sumu=sumu+u(i,T)
    return sumu/Z(T)

def U2(T):
    sumu2=0.0
    for i in range(len(E)):
        sumu2=sumu2+u2(i,T)
    return sumu2/Z(T)

def s(T):
    return np.sqrt(U2(T)-(U(T)*U(T)))
     
def cv(T):
     return s(T)**2/T**2

def F(T):
    return -kb*T*np.log(Z(T))

def S(T):
     return (U(T)-F(T))/T
def op(i,T):
   return np.exp(-E[i]/(kb*T))/Z(T)

par=[]
temp=[]
avg_e=[]
fluc=[]
sh=[]
free=[]
entropy=[]
Pr1=[]
Pr2=[]


for j in range(n):
    T=Ti+j*h
    temp.append(T)
    par.append(Z(T))
    avg_e.append(U(T))
    fluc.append(s(T))
    sh.append(cv(T))
    free.append(F(T))
    entropy.append(S(T))
    Pr1.append(op(0,T))
    Pr2.append(op(1,T))

plt.title("Internal Energy(U) vs Temperature(T) Plot",size=15)
plt.xlabel("Temperature($T$)", size=14)
plt.ylabel("Internal Energy($U$)",size=14)
plt.plot(temp,avg_e,marker='^',mfc='k',ms=3,ls='--',color='g',label='For Internal Energy($U$)')
plt.legend(loc='best',prop={'size':14}) 
plt.grid()
plt.show()
plt.title("Energy Fluctuation(ΔE) vs Temperature(T) Plot",size=15)
plt.xlabel("Temperature($T$)", size=14)
plt.ylabel("Energy Fluctuation ($ΔE$)",size=14)
plt.plot(temp,fluc,marker='*',mfc='k',ms=3,ls='--',color='purple',label='For Energy Fluctuation($ΔE$)')
plt.legend(loc='best',prop={'size':14}) 
plt.grid()
plt.show()
plt.title("Specific Heat(Cv) vs Temperature(T) Plot",size=15)
plt.xlabel("Temperature($T$)", size=14)
plt.ylabel("Specific Heat($C_v$)",size=14)
plt.plot(temp,sh,marker='.',mfc='r',ms=3,ls='--',color='r',label='For Specific Heat($C_v$)')
plt.legend(loc='best',prop={'size':14}) 
plt.grid()
plt.show()    
plt.title("Free Energy(F) vs Temperature(T) Plot",size=15)
plt.xlabel("Temperature($T$)", size=14)
plt.ylabel("Free Energy($F$)",size=14)
plt.plot(temp,free,marker='+',mfc='k',ms=3,ls='--',color='magenta',label='For Free Energy($F$)')
plt.legend(loc='best',prop={'size':14}) 
plt.grid()
plt.show()    
plt.title("Entropy(S) vs Temperature(T) Plot",size=15)
plt.xlabel("Temperature($T$)", size=14)
plt.ylabel("Entropy($S$)",size=14)
plt.plot(temp,entropy,marker='v',mfc='k',ms=3,ls='--',color='b',label='For Entropy($S$)')
plt.legend(loc='best',prop={'size':14}) 
plt.grid()
plt.show()
plt.title("Occupation number vs Temperature(T) Plot",size=15)
plt.xlabel("Temperature($T$)", size=14)
plt.ylabel("Occupation number($Pr_n$)",size=14)  
plt.plot(temp,Pr1,label='For occupation no.($Pr1$)')
plt.plot(temp,Pr2,label='For occupation no.($Pr2$)')
plt.legend(loc='best',prop={'size':14})
plt.grid()
plt.show()
