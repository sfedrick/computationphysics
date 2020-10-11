import numpy as np
import matplotlib.pyplot as plt
from pandas.core.common import flatten
from fourierMethods import *
from derivatives import *
def func1(x):
    y=np.sin(3*x)+3*np.cos(6*x)
    return y
def dx1(x):
    y=3*np.cos(3*x)-18*np.sin(6*x)
    return y
def func2(x):
    y=6*x-x**2
    return y
def dx2(x):
    y=6-2*x
    return y
N=32
x=np.linspace(0,2*np.pi,N)
K=findK(N)
y1=func1(x)
y2=func2(x)
A1=np.fft.fft(y1)/N
A2=np.fft.fft(y2)/N

'''
a1=A1[0:int(N/2)+1]

a2=A2[0:int(N/2)+1]


for i in range(1,len(a1)):
    a1[i]=a1[i]*K[i]*1j
    a2[i]=a2[i]*K[i]*1j
A1=realA(a1)
A2=realA(a2)
'''

for i in range(0,len(A1)):
    A1[i]=A1[i]*K[i]*1j
    A2[i]=A2[i]*K[i]*1j

print(len(A1))
Y1=np.fft.ifft(A1*N)
X1,D1=secondorderwhole(x,func1(x))
Y2=np.fft.ifft(A2*N)
X2,D2=secondorderwhole(x,func2(x))

error1spec=abs(Y1-dx1(x))
error1fd=abs(D1-dx1(X1))

error2spec=abs(Y2-dx2(x))
error2fd=abs(D2-dx2(X2))

figure1=plt.figure()
ax=figure1.add_subplot()
ax.plot(x,Y1,label="spectral",color='green')
ax.plot(X1,D1,label="second order FD",color='blue')
ax.plot(x,dx1(x),label="Exact",color='red')
ax.set_xlabel("x")
ax.set_ylabel("dy/dx")
ax.legend(loc='best')
figure1.suptitle("Function 1")

figure2=plt.figure()
ax=figure2.add_subplot()
ax.plot(x,Y2,label="spectral",color='green')
ax.plot(X2,D2,label="second order FD",color='blue')
ax.plot(x,dx2(x),label="Exact",color='red')
ax.set_xlabel("x")
ax.set_ylabel("dy/dx")
ax.legend(loc='best')
figure2.suptitle("Function 2")


figure3=plt.figure()
ax=figure3.add_subplot()
ax.plot(x,error1spec,label="spectral",color='green')
ax.plot(X1,error1fd,label="second order FD",color='blue')
ax.set_xlabel("x")
ax.set_ylabel("absolute error")
ax.legend(loc='best')
figure3.suptitle("error function 1")

figure4=plt.figure()
ax=figure4.add_subplot()
ax.plot(x,error2spec,label="spectral",color='green')
ax.plot(X2,error2fd,label="second order FD",color='blue')
ax.set_xlabel("x")
ax.set_ylabel("absolute error")
ax.legend(loc='best')
figure4.suptitle("error function 2")
plt.show()
