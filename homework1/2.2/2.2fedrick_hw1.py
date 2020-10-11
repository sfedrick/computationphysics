# -*- coding: utf-8 -*-
"""
Created on Sat Feb  1 15:28:22 2020

@author: shaun
"""
import numpy as np
#sets constants 
G=6.67*(10**(-11))
M=5.97*(10**(24))
R=6371.0
T=float(input("please enter the revolution time in minutes \n "))
T=T*60
#plugs in to derived formula
h=(((G*M*T**2)/(4*(np.pi)**2))**(1.0/3.0))-R*1000
print("the altitude needed for this revolution speed is: "+ str(h))