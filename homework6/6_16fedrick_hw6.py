from nonlinequationfedrick_hw6 import *
import matplotlib.pyplot as plt
import numpy as np
#function for secant method
def f(r):
    G=6.674*10**(-11)
    M=5.974*10**(24)
    m=7.348*10**(22)
    R=3.844*10**8
    w=2.662*10**(-6)
    try:
        y1=(G*M)/(r**2)
    except ZeroDivisionError:
        y1=0
    try:
        y2=(G*m)/((R-r)**2)
    except ZeroDivisionError:
        y2=0
    try:
        y3=r*w**2
    except ZeroDivisionError:
        y3=0
    y=y1-y2-y3
    return y

ans=secant(f,1e-6,1000,3.844*10**5)
print(ans)

R=np.linspace(1,3.844*10**8)
fig1=plt.figure()
ax=fig1.add_subplot()
ax.plot(R,f(R))
print(f(326045071.66535544))
plt.show()
