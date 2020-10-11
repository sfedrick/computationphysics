import numpy as np
import pandas as pd
from integration import *
filepathH = 'Hamiltonian.xls'
filepathVec='Eigenvectors.xlsx'
filepathVal='Eigenvalues.xlsx'
m=0
n=0
#length in meters
L=(5.0)*10.0**(-10)
#Length in angstroms
#L=5
# this was an epxiremental attempt to numerically solve the integral it doesn't work yet
def integrand(x):
    global m
    global n
    global L
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
    s=np.sin((np.pi*m*x)/L)
    v=a*x/L
    dx2=-(((np.pi*n*x)/L)**2)*np.sin((np.pi*n*x)/L)
    H=-((h**2)/(2*M))*dx2+v*np.sin((np.pi*n*x)/L)
    f=s*H
    return f


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
# size of matrix
N=10
# H matrix
Hamilt=np.empty([N,N])
'''
for row in range(0,N):
    for col in range(0,N):
        m=row+1
        n=col+1
        Hamilt[row][col]=swhole(100,0,L,integrand)[2]
'''
#fill H Matrix using the conditions found in part c
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
#get eigen values and eigen vectors of H matrix
evalues,evectors=np.linalg.eigh(Hamilt)
#get the ground state probability
ground=np.linalg.multi_dot([evectors[:,1],Hamilt,np.transpose(evectors[:,1])])
probsum=0
#sum over all states in the matrix
for i in range(1,N):
    probsum+=np.linalg.multi_dot([evectors[:,i],Hamilt,np.transpose(evectors[:,i])])
#outputs Matrices to an excel file for debugging
evalues=6.242*10**(18)*evalues
dfH = pd.DataFrame (Hamilt)
dfevec=pd.DataFrame (evectors)
dfevalues=pd.DataFrame(evalues)
dfH.to_excel(filepathH, index=False)
dfevec.to_excel(filepathVec, index=False)
dfevalues.to_excel(filepathVal, index=False)
