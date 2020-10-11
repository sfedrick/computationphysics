import numpy as np
import matplotlib.pyplot as plt
from pandas.core.common import flatten
from fourierMethods import *
from multipdesys import *
N=64
def g(x):
    y=2+6*np.sin(6*x)-38*np.cos(6*x)
    return y




x=np.linspace(0,2*np.pi,N)
G=g(x)
Fg=np.fft.fft(G)/N
FF=[]
K=findK(N)
for i in range(0,len(K)):
    ff=(-Fg[i])/(K[i]**2 + K[i]*1j +2)
    FF.append(ff)
FF=np.array(FF)*N
f=np.fft.ifft(FF)

d=np.copy(G)
D=d[1:-1]
L=len(D)
A=[]
B=[]
C=[]

dx=abs(x[0]-x[1])
alpha=(1/(dx**2))-(1/(dx))
beta=(-2/(dx**2))-2
gamma=(1/(dx**2))+(1/(dx))

for i in range(0,L):
    if i==0:
        C.append(alpha)
        B.append(beta)
    elif i==L-1:
        A.append(gamma)
        B.append(beta)
    else:
        A.append(alpha)
        B.append(beta)
        C.append(gamma)
xc=[0]
xc.append(list(TDMAsolver(A, B, C, D)))
xc.append([0])
xc=list(flatten(xc))



figure1=plt.figure()
ax=figure1.add_subplot()
ax.plot(x,f,label="spectral",color='red')
ax.plot(x,xc,label="FD",color='green')
ax.set_xlabel("x")
ax.set_ylabel("f(x)")
ax.legend(loc='best')
figure1.suptitle("Function 1")
plt.show()
