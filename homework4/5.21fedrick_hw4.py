# -*- coding: utf-8 -*-
"""
Created on Mon Feb 24 21:30:41 2020

@author: shaun
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import plotly.express as px
import plotly.graph_objects as go
def round_sig(x, sig):
    return round(x, sig-int(floor(log10(abs(x))))-1)

#a function to calculate the distance from a point in a grid
def distance(x,y,point,boxsize):
    #x component of r
    rx=(x-point[0])**2
    ry=(y-point[1])**2
    #y component
    r=rx+ry
    r=r**(0.5)
    #calculate r from the size of each bin in the grid
    r=r*boxsize
    return r
#make matrix it's odd so that the center is easy to find
N=101
#define the space to be from -0.5 meters to 0.5 meters in x and y
x=np.linspace(-0.5,0.5,N)
y=np.linspace(0.5,-0.5,N)
#create matrix A is the real matrix and Aplot is the matrix that is used to help in visualizing
A=np.empty([N,N],float)
Aplot=np.empty([N,N],float)
#center of the matrix
center=(N-1)/2
#in meters 
#sets the size of each grid in matrix A
boxsize=abs(x[0]-x[1])
#sets the location of the positive charge
a=np.array([center-(0.05/boxsize),center],float)
#sets the location of the negative charge
b=np.array([center+(0.05/boxsize),center],float)
#permitivity of freespace
e=8.85418782*(10**(-12))
#limit of acceptable potentials
limit=1e+010
#fill matrix with potetials at each point in the grid
for row in range(0,N):
    for col in range(0,N):
        potentiala=1/(4.0*np.pi*e*(distance(row,col,a,boxsize)))
        potentialb=-1/(4.0*np.pi*e*(distance(row,col,b,boxsize)))
        A[row][col]=(potentiala+potentialb)
        if((potentiala+potentialb)<-limit ):
            Aplot[row][col]=-limit
        elif((potentiala+potentialb)>limit):
            Aplot[row][col]=limit
        else:
            Aplot[row][col]=(potentiala+potentialb)
#plot the heat map
fig = go.Figure(data=go.Heatmap(
                   z=Aplot,
                   x=x,
                   y=y,
                   hoverongaps = False)
)

fig.update_layout(
    xaxis_title="X meters",
    yaxis_title="Y meters",
    title='Eletric Potential of two point charges'
    )
fig.show()
#find the gradient of A to get the eltric field plot using quiver
v, u = -np.gradient(Aplot, boxsize, boxsize)
A=np.log(abs(A))
figure1=plt.figure()
ax = figure1.add_subplot()
ax.set_xlabel("X meters")
ax.set_ylabel("Y meters")
q = ax.quiver(x, y, u, v)
figure1.suptitle("Eletric Field of 2 point charges")
plt.show()