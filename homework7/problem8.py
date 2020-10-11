from odesolversystem import *
import numpy as np
import matplotlib.pyplot as plt
from math import sin,pi
G=0
#defines the function system for part a
def af1(X,t):
    #ode system for the derivative of theta
    O1=X[0]
    O2=X[1]
    return O2
def af2(X,t):
    #outputs the second derivative of theta
    O1=X[0]
    O2=X[1]
    g=9.81
    l=0.1
    O=-g/l*sin(O1)
    return O
#these are the functions for driven system
def df1(X,t):
    O1=X[0]
    O2=X[1]
    return O2
def df2(X,t):
    O1=X[0]
    O2=X[1]
    g=9.81
    l=0.1
    C=2
    G=5
    O=-g/l*np.sin(O1)+ C*np.cos(O1)*sin(2*np.pi*G*t)
    return O
#this is the system for finding the resonant frequency
def Rf2(X,t):
    O1=X[0]
    O2=X[1]
    g=9.81
    l=0.1
    C=2
    O=-g/l*np.sin(O1)+ C*np.cos(O1)*sin(10*t)
    return O
#creates list of functions
aF=[af1,af2]
#creates list of init conditions
x0=np.array([179.0/180*np.pi,0])
#odesys ouputs a dictionary and a time vector
X,t=odesysRk4(x0,0.0,10,0.001,aF)
theta,omega=X.values()

#creates function list using driven functions
dF=[df1,df2]
#enters the init conditions
x0=np.array([0,0])
dX,dt=odesysRk4(x0,0.0,10.0,0.001,dF)
dtheta,domega=dX.values()


#creates the Gamma for your function this is a frequency
G=((9.8/0.1)**0.5)/(2*np.pi)
#creates function list using resonant function
rF=[df1,Rf2]
#enters the init conditions
x0=np.array([0,0])
rX,rt=odesysRk4(x0,0.0,100.0,0.001,rF)
rtheta,romega=rX.values()

#plots
figure=plt.figure()
ax=figure.add_subplot()
ax.plot(t,theta,label="Theta Position vs time")
ax.legend()
ax.set_xlabel("Time")
ax.set_ylabel("Theta Position radians")
figure.suptitle("Regular Pendulum rungekuta4")

figure1=plt.figure()
ax1=figure1.add_subplot()
ax1.plot(dt,dtheta)
ax1.set_xlabel("Time")
ax1.set_ylabel("Theta Position radians")
figure1.suptitle("Theta Position vs time Driven Pendulum C=2 Alpha=5")

figure2=plt.figure()
ax2=figure2.add_subplot()
ax2.plot(rt,rtheta)
ax2.set_xlabel("Time")
ax2.set_ylabel("Theta Position radians")
figure2.suptitle("Theta Position vs time Driven Pendulum C=2 Alpha="+str(G))

plt.show()
