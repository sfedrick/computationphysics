# -*- coding: utf-8 -*-
"""
Created on Tue Feb 18 15:11:00 2020

@author: shaun
"""

import numpy as np
from gaussxw import gaussxw
import matplotlib.pyplot as plt

N=50
#find legedre polynomial weights
x,w=gaussxw(N)
#define integrand function
def function(x):
    a=(x**4)*(np.e**(x))
    b=((np.e**(x))-1)**2
    y=a/b
    return y
#define heat capacity function
def f(t):
    d=6.022*(10**28)
    v=1000*(10**(-6))
    theta=428
    coeff1=(t/theta)**3
    coeff1=9*v*d*coeff1
    y=coeff1*Gintegral(0,theta/t,function)
    return y
#compute inegrand
def Gintegral(a,b,y):
    global N
    #rescale x and weights to the domain
    xp=0.5*(b-a)*x + 0.5*(b+a)
    wp=0.5*(b-a)*w
    s=0
    #find the value of the function at every x and multiply it by the weights to get the sum
    for xm in range(0,N):
        s+=wp[xm]*y(xp[xm])
    return s
#plot results
X=np.linspace(5,500)
y=f(X)
figure=plt.figure()
ax=figure.add_subplot()
ax.plot(X,y)
ax.set_xlabel("Temperature")
ax.set_ylabel(r"$C_{v}$")
figure.suptitle("Specific heat of Aluminum")







