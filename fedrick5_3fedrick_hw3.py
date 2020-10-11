# -*- coding: utf-8 -*-
"""
Created on Tue Feb 18 16:13:02 2020

@author: shaun
"""
import numpy as np
import matplotlib.pyplot as plt
from integration import *
N=1000
#creates function x
def function(x):
    y=np.e**(-(x**2))
    return y
#calculates the integral using n bins and simpsons rule
def f(x):
    global N
    return(swhole(N,0,x,function)[2])


X=np.linspace(0,3,30)
Y=[]
#plots result
for value in X:
    Y.append(f(value))
figure=plt.figure()
ax=figure.add_subplot()
ax.plot(X,Y)
ax.set_xlabel("X")
ax.set_ylabel(r"E(x)")
figure.suptitle("$E(x)$ evaluated numerically")
