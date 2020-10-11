# -*- coding: utf-8 -*-
"""
Created on Thu Jan 30 21:08:10 2020

@author: shaun
"""
import numpy as np
# Part A write a program that takes as it's inputs A and Z
A=float(input("enter value for A \n "))
Z=float(input("enter value for Z \n "))
# sets constants
a1=15.8
a2=18.3 
a3=0.714
a4=23.2

#checks init conditions to determine constants 
if (A%2!=0):
    a5=0
elif(A%2==0 and Z%2==0):
    a5=12.0
elif(A%2==0 and Z%2!=0):
    a5=-12.0
    

B=(a1*A-a2*A**(2.0/3.0)-a3*((Z**2)/(A**(1/3)))
   -a4*(((A-2*Z)**2)/A)+a5/(A**(0.5))
)
#outputs binding energy can check A=58 and Z=28 and it gives 
# a resonable result of 497.56 mv
print( "Given A="+str(A)+" Given Z="+str(Z)+" energy= " +str(B))
stable=0
stableA=0
#part B outbuts binding energy per nucleon
print( "binding energy per nucleon "+ str(B/A))
#part c 
#steps through z to 3z using 100 steps can decrease step size for more accuracy
listofenergies=np.linspace(Z,3*Z,100)

for x in listofenergies:
    A=x
    if (A%2!=0):
        a5=0
    elif(A%2==0 and Z%2==0):
        a5=12.0
    elif(A%2==0 and Z%2!=0):
        a5=-12.0
    B=(a1*A-a2*A**(2.0/3.0)-a3*((Z**2)/(A**(1/3)))
   -a4*(((A-2*Z)**2)/A)+a5/(A**(0.5))
   )
    if(stable<B/A):
        stable=B/A
        stableA=A
#prints out most stable value for inputted Z      
print("\n this is stable for given Z= "+ str(Z)+" \n"+ str(stable)+" the most\
      stable A is: "+ str(stableA))

#part d
stableglobal=0
stablez=0
epn=0
#runs through z=1 to z=100
for y in range(1,101):
    stable=0
    stableA=0
    Z=y
    #steps through every A between z and 3*z
    listofenergies=np.linspace(1,3*Z,num=100)
    for x in range(1,int(3*Z)):
        if(Z==28):
            print("this is the iteration "+str(x))
        A=x
        if (A%2!=0):
            a5=0
        elif(A%2==0 and Z%2==0):
            a5=12.0
        elif(A%2==0 and Z%2!=0):
            a5=-12.0
            
        B=(a1*A-a2*A**(2.0/3.0)-a3*((Z**2)/(A**(1/3)))
           -a4*(((A-2*Z)**2)/A)+a5/(A**(0.5))
           )
        if ((B/A)>stable):
            stable=B/A
            stableA=A
    if((B/A)>epn):
        stableglobal=B
        epn=B/A
        stablez=Z
    print("\n this is the stable A for Z= "+str(Z)+ " \n "+ str(stableA)\
          +" energy is equal to :"+ str(B)+ " B/A ="+ str(B/A))
            
            
print("\n this is greatest binding energy per nucleon "+ str(epn) +" and occurs at Z= "+ str(stablez))




