import numpy as np
from itertools import product

N,n=3,6
a=list(range(n+1))
l=list(product(a,repeat=N))

for j in l[:]:
    if sum (j) != n:
        l.remove(j)
        
k=np.array(l)        
print(k)    
print("no. of microstates are: ",len(l))    
'''from itertools import product

N,n=3,6
a=list(range(n+1))
l=list(product(a,repeat=N))

for j in l[:]:
    if sum (j) != n:
        l.remove(j)
        
print ("microstate no.        " ,"Energy 1  " , "Energy 2  " , "    Energy 3")
for i in range (len(l)):
   if i<9:
       print (i+1,"                    " , l[i][0],"             " , l[i][1],"             " , l[i][2])
   else :
       print (i+1,"                   " , l[i][0],"             " , l[i][1],"             " , l[i][2])'''
    