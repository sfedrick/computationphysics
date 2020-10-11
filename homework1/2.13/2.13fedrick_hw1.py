# -*- coding: utf-8 -*-
"""
Created on Sat Feb  1 17:37:06 2020

@author: shaun
"""
#Catalan numbers
def cn(n):
    if(n==0):
        return 1
    elif(n>0):
        return (((4*n-2)/(n+1))*cn(n-1))
#greatest common denominator 
def gcd (m,n):
    if(n==0):
        return m
    else:
        return gcd(n, m%n)

print("This is the 100th number in the catalan number series: "+str(cn(100)))

print("\n The greatest common divisor betwen 108 and 192: "+ str(gcd(192,108)))