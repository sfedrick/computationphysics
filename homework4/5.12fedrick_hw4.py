# -*- coding: utf-8 -*-
"""
Created on Mon Feb 24 17:48:21 2020

@author: shaun
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from sympy import *
from integration import *
from math import log10, floor


#redefine the basis so that you can integrate to infinity
def function(z):
    z=float(z)
    #calculate dz by taking the derivative
    dz=1/((1-z)**2)
    Z=(z/(1-z))
    #change xdx to be in the form zdz
    y=(Z**3)/((np.e**(Z))-1)*dz
    return y
#planks constant divided by 2pi
h=(6.62607004*(10**(-34)))/(2*np.pi)
#boltzman constant
k=1.38064852*(10**(-23))
#we set temperature to 1 so that it can be factored out 
T=1
#speed of light
c=299792458
#define the constant in front of the integral
bottom=4*(np.pi**2)*(c**2)*(h**3)
constant=((k**4)*(T**4))/bottom
x=[]
y=[]
#take the improper integral if you integrate to 0.999 you get a too large integer error
x,y,total=swhole(1000,0.000001,0.99,function)
#return stefan boltzman constant
ans=total*constant
print("This is the stefan boltzman constant "+ str(ans))