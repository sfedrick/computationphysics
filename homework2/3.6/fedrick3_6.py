# -*- coding: utf-8 -*-
"""
Created on Tue Feb 11 08:46:36 2020

@author: shaun
"""
import numpy as np
import matplotlib as mp
#creates recursive function for logisitic map 
def logmap(n,r,x):
    if(n==0):
        return r*x*(1-x)
    else:
        return logmap(n-1,r,r*x*(1-x))
#create linear space from 1 to 4 with step size 0.01
r=np.linspace(1,4,300)
y=[]
y2=[]
y3=[]
y4=[]   
#fill in y values 
for i in r:
    #in order to capture the oscllatory motion you have to change the iterations 
    #to be 1000, 10001 10002 and 10003
    y.append(logmap(1000,i,0.5))
    y2.append(logmap(1001,i,0.5))
    y3.append(logmap(1002,i,0.5))
    y4.append(logmap(1003,i,0.5))

#create plot 
fig1=mp.pyplot.figure()
ax=fig1.add_subplot(111)
fig1.suptitle("Figtree")
ax.scatter(r,y,c='k', marker="o", label='Original Data')
ax.scatter(r,y2,c='k', marker="o", label='Original Data')
ax.scatter(r,y3,c='k', marker="o", label='Original Data')
ax.scatter(r,y4,c='k', marker="o", label='Original Data')
ax.set_xlabel(r"R value $\Delta$ 0.01")
ax.set_ylabel("F(x)")


