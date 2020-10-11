# -*- coding: utf-8 -*-
"""
Created on Tue Feb 11 00:14:17 2020

@author: shaun
"""
import matplotlib as mp 
import numpy as np
#create the deltoid function using catesian coordinates
def deltoid(theta):
    x=2*np.cos(theta)+np.cos(2*theta)
    y=2*np.sin(theta)-np.sin(2*theta)
    return x,y
#create deltoid function using polar coordinates
def deltoidpolar(theta):
    r=theta**2
    x=r*np.cos(theta)
    y=r*np.sin(theta)
    return x,y
#create the fey function using polar coordinates
def feyfunction(theta):
    r=((np.e)**(np.cos(theta)))-2*np.cos(4*theta)+ \
        np.sin(theta/12)**5
    x=r*np.cos(theta)
    y=r*np.sin(theta)
    return x,y

#initialize list for deltoid function
xvals=[]
yvals=[]
step=np.linspace(0,2*np.pi, 100)
# fill in list for deltoid function
for i in step:
    x,y=deltoid(i)
    xvals.append(x)
    yvals.append(y)
    
fig1=mp.pyplot.figure()
ax=fig1.add_subplot(1,1,1)
ax.plot(xvals, yvals, c='b', marker="s", label='Original Data')
fig1.suptitle("Deltoid")
ax.set_xlabel(r"Angle from 0 to  $2\pi$")
ax.set_ylabel("F(x)")

#initialize polar deltoid lists
xpol=[]
ypol=[]
step2=np.linspace(0,10*np.pi, 100)
# fill in list for polar deltoid function 
for i in step2:
    x,y=deltoid(i)
    xpol.append(x)
    ypol.append(y)
 #plot   
fig2=mp.pyplot.figure()
ax1=fig2.add_subplot(1,1,1)
ax1.plot(xpol, ypol, c='b', marker="s", label='Original Data')
fig2.suptitle("Deltoid Polar form")
ax1.set_xlabel(r"Angle from 0 to $10\pi$")
ax1.set_ylabel("F(x)")
 
    
    
#initialize fey lists   
xfey=[]
yfey=[]   
step3=np.linspace(0,24*np.pi, 5000 )
#fill in fey lists
for i in step3:
    x,y=feyfunction(i)
    xfey.append(x)
    yfey.append(y)
#plot fey function
fig3=mp.pyplot.figure()
ax2=fig3.add_subplot(1,1,1)
ax2.plot(xfey, yfey, c='b', marker="s", label='Original Data')
fig3.suptitle("Fey's Function")
ax2.set_xlabel("Angle from 0 to $24\pi$")
ax2.set_ylabel("F(x)")
     
    
    
    
    
    
    
    
    
    
    
    
    
    
    