# -*- coding: utf-8 -*-
"""
Created on Sun Feb 16 13:25:30 2020

@author: shaun
"""

import numpy as np
import matplotlib.pyplot as plt
from sympy import *

#everything with part computes the integral of a single bin the letter after signifies
#which method it uses
#r is rectangular
#t is trapezoidal
#s is simpsons
def rpart(y,delta):
    sumx=y*delta
    return sumx

def tpart(y1,y2,delta):
    sumx=((y1+y2)/2.0)*delta
    return sumx

def spart(y1,y2,y3,delta):
    sumx=(float(delta)/3.0)*(y1+(4.0*y2)+y3)
    return sumx

#everything with whole adds al the bins together and computes the total integral from a to b using
#N points like with part the letter signifies the type of method it uses

def rwhole(N,a,b,left,function):
    listx=np.linspace(a,b,N)
    listy=[]
    sumx=0
    if(left==True):
        for x in range(0,N):
            if(x<N-1):
                y1=function(listx[x])
                delta=listx[x]-listx[x+1]
                delta=(delta**2)**(0.5)
                sumx+= rpart(y1,delta)
                listy.append(sumx)
            else :
                return sumx
    else:
         for x in range(0,N):
            if(x<N-2):
                y2=function(listx[x+1])
                delta=listx[x]-listx[x+1]
                delta=(delta**2)**(0.5)
                sumx+= rpart(y2,delta)
                listy.append(sumx)
            else :
                return listx,listy,sumx
def twhole(N,a,b,function):
    listx=np.linspace(a,b,N)
    listy=[]
    newlistx=[]
    sumx=0
    for x in range(0,N-1):
        y1=function(listx[x])
        y2=function(listx[x+1])
        delta=listx[x]-listx[x+1]
        delta=(delta**2)**(0.5)
        sumx+= tpart(y1,y2,delta)
        newlistx.append(listx[x])
        listy.append(sumx)
    return newlistx,listy,sumx

#this evaluates the trapezoidal rule with end corrections
def twhole2(N,a,b,function):
    listx=np.linspace(a,b,N)
    delta=listx[0]-listx[1]
    listy=[]
    sumx=0
    A=fderiv(a)
    B=fderiv(b)
    sumx-=((delta**2)/12.0)*(B-A)
    for x in range(0,N):
        if(x<N-1):
            y1=function(listx[x])
            y2=function(listx[x+1])
            delta=(delta**2)**(0.5)
            sumx+= tpart(y1,y2,delta)
            listy.append(sumx)
        else :
            return listx,listy,sumx


def swhole(N,a,b,function):
    listx=np.linspace(a,b,N)
    delta=listx[0]-listx[1]
    listy=[]
    sumx=0
    for x in range(1,int(N/2)):
        y1=function(listx[2*x-2])
        y2=function(listx[2*x-1])
        y3=function(listx[2*x])
        delta=(delta**2.0)**(0.5)
        sumx+= spart(y1,y2,y3,delta)
        listy.append(sumx)
    return listx,listy,sumx
# rhomberg integration implements smaller and smaller step sizes to take a linear combination of an
# integral technique with a bin size of h and several other bin sizes of h/n
def rhomberg(N,a,b,steps):
    L=2**N
    h=(float(b)-(a))/N
    A=zeros(2**N)
    B=zeros(2**N,1)
    print("We are currently at bin= "+ str(b))
    for x in range(0,L):
        for y in range(0,L):
            if y==0:
                A[x,y]=1.0
                n=(2**x)*steps
                print("row= "+str(x) +" has this many steps "+ str(n))
                B[x,0]=swhole(int(n),a,b)[2]

            else:
                A[x,y]=-(h/(2**x))**(2*y)
    #computes the combination of integral bins to reduce the error of the integral
    System=(A.copy()).col_insert((L),B)
    x=symbols('a0:%d'%(L),Real=True)
    b=solve_linear_system(System,*x)

    evals=0
    i=N
    while i>0:
        evals= 2**(i)+evals
        i=i-1
    ans=[b[x[0]],b[x[L-1]],evals]

    return ans

#uses rhomberg part at each location to compute the entirety of the rhomberg integral
def rhombergwhole(a,b,step,error):
    bins=np.linspace(a,b,step)
    x=[]
    y=[]
    evals=[]
    sumy=0;
    for i in range(0,len(bins)-1):
        exitwhile=False
        locala=bins[i]
        localb=bins[i+1]
        N=1
        while( not exitwhile):
            ans=rhomberg(N,locala,localb,step)
            #print("this is the N value "+str(N))
            #print("this is the error " + str(ans[1]))
            if(ans[1]<error):
                exitwhile=True
            if(N>2):
                exitwhile=True
            N=N+1
        evals.append(ans[2])
        sumy=sumy+float(ans[0])
        x.append(bins[i])
        y.append(sumy)
    return x,y,evals
