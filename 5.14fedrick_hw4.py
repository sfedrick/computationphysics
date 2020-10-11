# -*- coding: utf-8 -*-
"""
Created on Mon Feb 24 18:39:17 2020

@author: shaun
"""

import numpy as np
from gaussxw import gaussxw
import matplotlib.pyplot as plt
N=100
#grab key points and weights from legendre polymials
x,w=gaussxw(N)
#define the integrand with input x,y,z
def integrand(x,y,z):
    f=1/(((x**2)+(y**2)+(z**2))**(3.0/2))
    return f
#define the function z 
def function(z):
    density=(10.0*1000)/(100)
    G=6.674*(10**(-11))
    F=G*density*z*Gmulti(-5,5,-5,5,integrand,z)
    return F
#take the triple integral
def Gmulti(a,b,c,d,f,z):
    global N
    global x
    global w
    #rescale x and weights to the domain
    xp=0.5*(b-a)*x + 0.5*(b+a)
    wpx=0.5*(b-a)*w
    yp=0.5*(d-c)*x + 0.5*(d+c)
    wpy=0.5*(d-c)*w
    s=0
    
    for ym in range(0,N):
        #find the value of the function at every y and multiply it by the weights to get the sum
        for xm in range(0,N):
            #find the value of the function at every x and multiply it by the weights to get the sum
            s+=wpy[ym]*wpx[xm]*f(xp[xm],yp[ym],z)
    return s


   

#plot results
Z=np.linspace(0,10,100)
F=function(Z)
figure=plt.figure()
ax=figure.add_subplot()
ax.plot(Z,F)
ax.set_xlabel("Z")
ax.set_ylabel(r"$Force_{z}$")
figure.suptitle("Force in Z direction")

plt.show()
      