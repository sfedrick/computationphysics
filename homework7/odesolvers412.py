# -*- coding: utf-8 -*-
"""
Created on Wed Feb 26 19:24:02 2020

@author: shaun
"""
import numpy as np
from scipy.optimize import brentq
def eulerstep(yn,tn,f,h):
    yn1=yn+h*f(yn,tn)
    return yn1
def eIlinstep(yn,tn,f,df,h):
    top=yn+h*f(yn,tn+h)-h*yn*df(yn,tn+h)
    bot=1-h*df(yn,tn+h)
    yn1=top/bot
    return yn1
def nonlinE(yn,tn1,h,f):
	#solve the differential equation you have for
	# the equation that satisfies yn+1 equaling some function of x and t
    c=yn
    def function(x):
        return -x+c+np.e**(x-tn1)
    n=np.linspace(-100,100,100)
    a=0
    b=0
    for elem in range(0,len(n)-1):
        A=function(n[elem])
        B=function(n[elem+1])
        if(A<0 and B>0):
            a=n[elem]
            b=n[elem+1]
        elif(A>0 and B<0):
            a=n[elem]
            b=n[elem+1]
    try:
        root=brentq(function,a,b)
    except:
        print("an error occured  will use euler for value \n")
        print("c= "+str(c))
        print("tn1= "+str(tn1))
        root=eulerstep(yn,tn1,f,h)
    return root
def trapezoidalE(yn,tn1,tn,h):
	#solve the differential equation you have for yn+1
	A= yn+(h/2.0)*2*((1+tn1)**3)*(np.e**(-tn1))
	B=(h/2.0)*(-3*tn*yn)/(1+tn1)
	C=(h/2.0)*2*((1+tn)**3)*(np.e**(-tn))
	top=A+B+C
	bottom=1+(h/2.0)*((3*tn1)/(1+tn1))
	yn1=top/bottom
	return yn1
def rk2step(yn,tn,f,h):
    y1=yn +(h/2)*f(yn,tn)
    yn1=yn+h*(y1,tn+tn*1/2)
    return yn1
def rk4step(yn,tn,f,h):
    a=1/2.0
    k1=h*f(yn,tn)
    k2=h*f(yn+a*k1,tn+h/2.0)
    k3=h*f(yn+a*k2,tn+h/2.0)
    k4=h*f(yn+a*k3,tn+h)
    yn1=yn+(1/6)*k1+(1/3)*(k2+k3)+(1/6)*k4

def eulerE(y0,a,b,f,h):
    y=[y0]
    t=[a]
    start=a
    end=a
    while(end<b):
        newy=eulerstep(y[-1],t[-1],f,h)
        y.append(newy)
        t.append(t[-1]+h)
        end+=h
    return t,y

def eulerI(y0,a,b,h,f):
    y=[y0]
    t=[a,a+h]
    start=a
    end=a
    while(end<b):
        newy=nonlinE(y[-1],t[-1],h,f)
        y.append(newy)
        t.append(t[-1]+h)
        end+=h
    del t[-1]
    return t,y
def eulerIlin(y0,a,b,h,f,df):
    y=[y0]
    t=[a]
    start=a
    end=a
    while(end<b):
        newy= eIlinstep(y[-1],t[-1],f,df,h)
        y.append(newy)
        t.append(t[-1]+h)
        end+=h
    return t,y
def eulerT(y0,a,b,h):
    y=[y0]
    t=[a,a+h]
    start=a
    end=a
    while(end<b):
        newy=trapezoidalE(y[-1],t[-1],t[-2],h)
        y.append(newy)
        t.append(t[-1]+h)
        end+=h
    del t[-1]
    return t,y
def rk2(y0,a,b,f):
    y=[y0]
    t=[a]
    start=a
    end=a
    while(end<b):
        newy=rk2step(y[-1],t[-1],f,h)
        y.append(newy)
        t.append(t[-1]+h)
        end+=h
    return t,y
def rk4(y0,a,b,h,f):
    y=[y0]
    t=[a]
    start=a
    end=a
    while(end<b):
        newy=rk4step(y[-1],t[-1],f,h)
        y.append(newy)
        t.append(t[-1]+h)
        end+=h
    return t,y
