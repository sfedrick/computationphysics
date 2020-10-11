# -*- coding: utf-8 -*-
"""
Created on Tue Feb 18 10:28:58 2020

@author: shaun
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from sympy import *
from integration import *
from math import log10, floor

#set global variables that are changed within the function bessel 
# this avoids having to pass more parameters into bessel
besselx=0
M=0
#define integrand within bessel function
def function(x):
    global M
    global besselx
    global test
    function=np.cos(M*x-besselx*np.sin(x))  
    test=1000
    return function
#define bessel function
def bessel(m,x):
    global besselx
    global M
    besselx=x
    M=m
    f=(1/np.pi)*swhole(1000,0,np.pi,function)[2]
    return f
#rounds to the appropriate sig fig when calculating density
def round_sig(x, sig):
    return round(x, sig-int(floor(log10(abs(x))))-1)

listx=[]
listy=[]
#creates a list of bessel functions with different m values
for m in range(0,20):
    listy.append([])
    listx.append(m)
    for x in range(0,20):
        listy[-1].append(bessel(m,x))
#plot these bessel functions with different symbols and colors 
figure1=plt.figure(figsize=(20,10))
for y in range(0, len(listy)):
    color="C"+str(y)
    m=y%12
    plt.plot(listx,listy[y],color,marker=m, label="M="+str(y))
#plot
plt.legend()    
plt.xlabel("X")
plt.ylabel("Value of bessel function")
plt.title("Bessel function at different values of M")



#creat array for heat map

A=np.zeros([100,100],float)
row=A.shape[0]
col=A.shape[1]
xspace=np.linspace(-1,1,row)
yspace=np.linspace(-1,1,col)
intensitymax=0
#walk through array and sets density at each point in array
for x in range(0,row):
    print("We are "+str(100*float(x)/row)+" percent finished")
    for y in range(0,col):
        xvalue=xspace[x]
        yvalue=yspace[y]
        #calculate r
        r=(xvalue**2+yvalue**2)**0.5
        lamb=0.5
        k=(2*np.pi)/lamb
        if(r!=0):
            j=bessel(1,k*r)
            intensity=(j/(k*r))**2
            if (intensity>intensitymax):
                intensitymax=intensity
        else:
            intensity=1/4.0
        A[x][y]=round_sig(intensity,1)
#plot the heat map
figure2=plt.figure(figsize=(5,5))
ax=figure2.add_subplot(1,1,1)
im=ax.imshow(A,vmax=intensitymax)      
figure2.suptitle(r"Light intensity Density plot",fontsize=30)
ax.set_xlabel(r"X",fontsize=15)
ax.set_ylabel(r"Y",fontsize=15)
# create legend for heat map
v = np.unique(A.ravel())
#include 0 in the heat map legend
values=[0]
for value in v:
    #only shows values greater than 0.01 in the heat map legend 
    if value>0.01:
       values.append(value)  
colors = [im.cmap(im.norm(value)) for value in values] 
patches = [ matplotlib.patches.Patch(color=colors[i], label="Density {l}".format(l=values[i]) ) for i in range(len(values)) ]
ax.legend(handles=patches, bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0. )
plt.show()   
    
