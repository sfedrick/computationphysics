import numpy as np
import pandas as pd
from integration import *
import matplotlib.pyplot as plt
filepathH = 'Hamiltoniansimple.xlsx'
filepathVec='Eigenvectorssimple.xlsx'
filepathVal='Eigenvaluessimple.xlsx'
m=0
n=0
#length in meters
L=(5.0)*10.0**(-10)
#Length in angstroms
#L=5


#mass in kg
M=9.1094*10**-31
#charge in joules
a=10*1.602176634*10**(-19)
#charge in eV
#a=10
#hbar joules per second
h=1.054571817*10**(-34)
#hbar eV per second
#h=6.582119569*10**(-16)

N=10
Hamilt=np.empty([N,N])

width=np.linspace(0,(5.0)*10.0**(-10),100)
e1=[]
e2=[]
e3=[]
#walks through the width of the box and adds the probability at x in the box
for L in width:
    for row in range(0,N):
        for col in range(0,N):
            m=row+1
            n=col+1
            if(m==n):
                Hamilt[row][col]=(((M*a*L**2)+((n**2)*(np.pi**2)*(h**2))))/(2*M*L**2)
            elif((m%2!=0 and n%2==0) or (m%2==0 and n%2!=0)):
                Hamilt[row][col]=((-8*a)/(np.pi**2))*((m*n)/(m**2 -n**2)**2)
            elif(m%2==0 and n%2==0) or (m%2!=0 and n%2!=0):
                Hamilt[row][col]=0

    evalues,evectors=np.linalg.eigh(Hamilt)
    e1.append(np.linalg.multi_dot([evectors[:,1],Hamilt,np.transpose(evectors[:,1])]))
    e2.append(np.linalg.multi_dot([evectors[:,2],Hamilt,np.transpose(evectors[:,2])]))
    e3.append(np.linalg.multi_dot([evectors[:,3],Hamilt,np.transpose(evectors[:,3])]))
#plot figure
figure1=plt.figure()
a1=figure1.add_subplot()
a1.plot(width,e1,color='r',label='Ground state')
a1.plot(width,e2,color='b',label='First excited state')
a1.plot(width,e3,color='g',label='Second excited state')
figure1.suptitle("Probability distributions for each state")
a1.set_ylabel("Probability")
a1.set_xlabel("L")
a1.legend()
plt.show()
