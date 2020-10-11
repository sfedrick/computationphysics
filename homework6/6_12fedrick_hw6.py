import pandas as pd
import  scipy as sp
import matplotlib.pyplot as plt
from io import StringIO
import numpy as np
from scipy.interpolate import interp1d
from derivativesfedrick_hw6 import *
from integrationfedrick_hw6 import *
from nonlinequationfedrick_hw6 import *
#function for relaxation method rearranged
def f(X):
    a=1
    b=2
    x=X[0]
    y=X[1]
    A=[]
    for l in X:
        A.append(0)
    A[0]=((b/y)-a)**0.5
    A[1]=x/(a + x**2)
    #X[1]=((a+b**2)/b)*(y**2)
#function for newton raphson method
def f1(X):
    a=1
    b=2
    x=X[0]
    y=X[1]
    A=[]
    for l in X:
        A.append(0)
    A[0]=((b/y)-a)**0.5-x
    A[1]=x/(a + x**2)-y
    #X[1]=((a+b**2)/b)*(y**2)
    return A
#computes the jacobian matrix for the newton raphson system
def J(X):
    a=16
    b=2
    x=X[0]
    y=X[1]
    f1dxa=2*(x**3)*y+2*x*y*(a+x**2)
    f1dxb=(y*(a+x**2))**2
    #f1dx=(f1dxa/f1dxb)-1
    f1dx=-1
    #f1dy=(-(x**2))/((y**2)*(a+x**2))
    f1dy=0
    row1=[f1dx,f1dy]
    #f2dx=(2*x*(y**2))/b
    f2dx=((a+x**2)-(2*x**2))/((a+x**2)**2)
    #f2dy=2*((a+x**2)/b)*y-1
    f2dy=-1
    row2=[f2dx,f2dy]
    j=np.array([row1,row2])
    return j
#output interation and results
r,ri=relax([-100,100],f, 1e-6)
nr, nri=NRsys(f1,J,1e-6,[-100,100])

print("this is the solution using relaxation: "+ str(np.array(r).real)+" with "+str(ri)+" iterations")
print("this is the solution using newton raphson: "+ str(np.array(nr).real)+" with "+str(nri)+" iterations")
