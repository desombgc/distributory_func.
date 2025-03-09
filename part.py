import numpy as np
import matplotlib.pyplot as plt

Ti,Tf=0.01,20.00                 #temperature range
n=100                         #no. of points
h=(Tf-Ti)/n                   
kb=1.0                        #boltzman constant taken as unity.
e1=1.0                        
e2=2.0

#defining Partition Function
def Z(T):
    return np.exp(-e1/(kb*T))+np.exp(-e2/(kb*T))
#defining Average Energy
def U(T):
    return ((e1*np.exp(-e1/(kb*T)))+(e2*np.exp(-e2/(kb*T))))/Z(T)
#defing Average Square Energy
def U2(T):
    return (((e1*e1)*np.exp(-e1/(kb*T)))+((e2*e2)*np.exp(-e2/(kb*T))))/Z(T)
#defining Energy fluctuation 
def s(T):
    return np.sqrt(U2(T)-(U(T)*U(T)))
#defing Specific Heat at constant volume     
def cv(T):
     return s(T)**2/T**2
#definig Free Energy
def F(T):
    return -kb*T*np.log(Z(T))
#defining Entropy
def S(T):
     return (U(T)-F(T))/T
#defining Ocuupation Number of two systems
def op1(T):
    return np.exp(-e1/(kb*T))/Z(T)

def op2(T):
    return np.exp(-e2/(kb*T))/Z(T)

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
    Pr1.append(op1(T))
    Pr2.append(op2(T))
##Plotting partition function(Z) vs Temperature(T)
plt.title("Partition Function(Z) vs Temperature(T) Plot",size=15)
plt.xlabel("Temperature($T$)", size=14)
plt.ylabel("Partition Function($Z$)",size=14)
plt.plot(temp,par,marker='o',mfc='k',ms=3,ls='--',color='b',label='For Partition Function($Z$)')
plt.legend(loc='best',prop={'size':14})
plt.show()

##1.A)Plotting internal energy(U) as a function of kbT
plt.title("Internal Energy(U) vs Temperature(T) Plot",size=15)
plt.xlabel("Temperature($T$)", size=14)
plt.ylabel("Internal Energy($U$)",size=14)
plt.plot(temp,avg_e,marker='^',mfc='k',ms=3,ls='--',color='g',label='For Internal Energy($U$)')
plt.legend(loc='best',prop={'size':14}) 
plt.show()

#1.B)Plotting Energy Fluctuation as a function of kbT
plt.title("Energy Fluctuation(ΔE) vs Temperature(T) Plot",size=15)
plt.xlabel("Temperature($T$)", size=14)
plt.ylabel("Energy Fluctuation ($ΔE$)",size=14)
plt.plot(temp,fluc,marker='*',mfc='k',ms=3,ls='--',color='purple',label='For Energy Fluctuation($ΔE$)')
plt.legend(loc='best',prop={'size':14}) 
plt.show()

##1.C)Plotting Specific Heat (Cv) at constant volume as a function of kbT
plt.title("Specific Heat(Cv) vs Temperature(T) Plot",size=15)
plt.xlabel("Temperature($T$)", size=14)
plt.ylabel("Specific Heat($C_v$)",size=14)
plt.plot(temp,sh,marker='.',mfc='r',ms=3,ls='--',color='r',label='For Specific Heat($C_v$)')
plt.legend(loc='best',prop={'size':14}) 
plt.show()

##1.D)Plotting Free Energy(F) as a function of kbT
plt.title("Free Energy(F) vs Temperature(T) Plot",size=15)
plt.xlabel("Temperature($T$)", size=14)
plt.ylabel("Free Energy($F$)",size=14)
plt.plot(temp,free,marker='+',mfc='k',ms=3,ls='--',color='magenta',label='For Free Energy($F$)')
plt.legend(loc='best',prop={'size':14}) 
plt.show()

##1.D)Plotting Entropy(S) as a function of kbT
plt.title("Entropy(S) vs Temperature(T) Plot",size=15)
plt.xlabel("Temperature($T$)", size=14)
plt.ylabel("Entropy($S$)",size=14)
plt.plot(temp,entropy,marker='v',mfc='k',ms=3,ls='--',color='b',label='For Entropy($S$)')
plt.legend(loc='best',prop={'size':14}) 
plt.show()

##Plotting Occupation number(op) as a function of kbT
plt.title("Occupation number vs Temperature(T) Plot",size=15)
plt.xlabel("Temperature($T$)", size=14)
plt.ylabel("Occupation number($Pr_n$)",size=14)
plt.plot(temp,Pr1,marker='o',mfc='k',ms=3,ls='-.',color='r',label='For occupation no.($Pr1$)')
plt.plot(temp,Pr2,marker='o',mfc='g',ms=3,ls='-.',color='b',label='For occupation no.($Pr2$)')
plt.legend(loc='best',prop={'size':14}) 
plt.show()