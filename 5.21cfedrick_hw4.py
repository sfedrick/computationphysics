# -*- coding: utf-8 -*-
"""
Created on Tue Feb 25 00:16:07 2020

@author: shaun
"""



import numpy as np
from gaussxw import gaussxw
import matplotlib.pyplot as plt
import plotly.express as px
import plotly.graph_objects as go
#grab key points and weights from legendre polymials
N=10
x,w=gaussxw(N)
#define the integrand with insput x,y and x0 and y0
#x0 and y0 indicate your location in the matrix
def integrand(x,y,x0,y0):
    q=100
    density=q*np.sin(2*np.pi*x/0.05)*np.sin(2*np.pi*y/0.05)
    f=density/(((x-x0)**2+(y-y0)**2)**0.5)
    return f

def function(x0,y0):
    #sets the bound of charge distribution to be a 10 cm by 10 cm square in the center of the matrix
    F=Gmulti(-0.05,0.05,-0.05,0.05,integrand,x0,y0)
    return F
#take the triple integral
def Gmulti(a,b,c,d,f,x0,y0):
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
            s+=wpy[ym]*wpx[xm]*f(xp[xm],yp[ym],x0,y0)
    return s
#set size of matrix 
size=100
#define the range and domain
X=np.linspace(-0.5,0.5,size)
Y=np.linspace(0.5,-0.5,size)
A=np.empty([size,size],float)
#specify the size of matrix grids
boxsize=abs(X[1]-X[0])
limit=0.05
#fill A with the potential given the location of the grid 
for row in range(0,len(Y)):
    print("percent done "+ str(float(row)/len(Y)))
    for col in range(0,len(X)):
        if(-limit<X[col]<limit and -limit<Y[row]<limit):
             A[row][col]=0
        else:
            A[row][col]=function(X[col],Y[row])

#create heat map
fig = go.Figure(data=go.Heatmap(
                   z=A,
                   x=X,
                   y=Y,
                   zmin=-0.01, 
                   zmax=0.01,
                   zauto=False,
                   hoverongaps = False)                
)

fig.update_layout(
    xaxis_title="X meters",
    yaxis_title="Y meters",
    title='Eletric Potential of continous charge distribution'
    )
fig.show()
#find the gradient of A to get the eltric field plot using quiver
v, u = np.gradient(A, boxsize, boxsize)
figure1=plt.figure()
ax = figure1.add_subplot()
ax.set_xlabel("X meters")
ax.set_ylabel("Y meters")
q = ax.quiver(X, Y, u, v)
figure1.suptitle("Eletric Field of continous charge distributions")
plt.show()