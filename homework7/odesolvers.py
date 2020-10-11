# -*- coding: utf-8 -*-
"""
Created on Wed Feb 26 19:24:02 2020

@author: shaun
"""
import numpy as np
def eulerstep(yn,tn,f,h):
    yn1=yn+h*f(yn,tn)
    return yn1
def nonlinE(yn,tn1,h):
	#solve the differential equation you have for
	# the equation that satisfies yn+1 equaling some function of x and t
	top=yn+h*2.0*((1.0+tn1)**3.0)*np.e**(-tn1)
	bottom=1.0+h*((3.0*tn1)/(1.0+tn1))
	yn1=top/bottom
	return yn1
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
    y1=yn+(h/2)*f(yn,tn)
    yn1=yn+h*f(y1,tn+tn*1/2)
    return yn1
def rk4step(yn,tn,f,h):
    a=1/2.0
    k1=h*f(yn,tn)
    k2=h*f(yn+a*k1,tn+h/2.0)
    k3=h*f(yn+a*k2,tn+h/2.0)
    k4=h*f(yn+a*k3,tn+h)
    yn1=yn+(1/6)*k1+(1/3)*(k2+k3)+(1/6)*k4
    return yn1
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

def eulerI(y0,a,b,h):
    y=[y0]
    t=[a,a+h]
    start=a
    end=a
    while(end<b):
        newy=nonlinE(y[-1],t[-1],h)
        y.append(newy)
        t.append(t[-1]+h)
        end+=h
    del t[-1]
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
def rk2(y0,a,b,f,h):
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
def rk4(y0,a,b,f,h):
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
