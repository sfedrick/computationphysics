from odesolversystem import *
import numpy as np
import matplotlib.pyplot as plt
from math import sin,pi
#create functions for the ode
def f1(X,t):
    x=X[0]
    y=X[1]
    dx=x-0.5*x*y
    return dx
def f2(X,t):
    x=X[0]
    y=X[1]
    dy=0.5*x*y-2*y
    return dy
#set init conditions
x0=[2,2]
F=[f1,f2]
#input values into the odesys
X,t=odesysRk4(x0,0,30,0.001,F)
#unpack values from the dictionary
rabbits,foxes=X.values()
figure=plt.figure()
#plot these valu
ax=figure.add_subplot()
ax.plot(t,rabbits,label="rabbits")
ax.plot(t,foxes,label="foxes")
ax.legend()
ax.set_xlabel("Time")
ax.set_ylabel("Population(Thousands)")
figure.suptitle("Predators and Prey")
plt.show()
