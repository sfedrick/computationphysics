import numpy as np
import matplotlib.pyplot as plt
from pandas.core.common import flatten
from fourierMethods import *
def function(x):
    y=np.cos(2.0*x)+0.5*np.cos(4.0*x)+(1.0/6)*np.cos(12.0*x)
    return y
T=2*np.pi
N=8
x=np.linspace(0,2*np.pi,N)
y=function(x)
k=np.fft.fftfreq(x.shape[-1])
K=findK(N)
print(K)
A=np.fft.fft(y)/N

NewA=abs(A)
B=np.fft.ifft(N*A)

figure1=plt.figure()
ax=figure1.add_subplot()
ax.scatter(K,NewA)
ax.set_xlabel("K")
ax.set_ylabel("A(k)")

figure2=plt.figure()
ax2=figure2.add_subplot()
ax2.plot(x,y,label="Original Function")
ax2.scatter(x,B,label="Inverse transformed function")
ax2.legend(loc='best')
ax2.set_xlabel("X")
ax2.set_ylabel("F(x)")

plt.show()
