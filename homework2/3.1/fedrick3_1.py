# -*- coding: utf-8 -*-
"""
Created on Mon Feb 10 13:15:55 2020

@author: shaun
"""
from numpy import loadtxt as lt
from matplotlib import pyplot as plt
#crates run average function to compute the running
#average given a list of 11 points 
def runaverage(list5):
    r=5
    a= 1/(2*r+1)
    sum5=0
    for x in list5:
        sum5+=x
    i=a*sum5
    return i
 
#read in data 
data=lt(r"C:\Users\shaun\Programs\Python\Computational Physics\cpresources\sunspots.txt",float)
ox=data[:,0]
y=data[:,1]
# make sure the inputed number of points is an integer
try:
    points=int(input("please enter the number of points you'd like to plot \n"))
except ValueError:
    print("you didn't enter a valid number try again next time we'll set the number of points to 1000")
    points=1000
#select the given number of points
newx=ox[0:points]
newy=y[0:points]
running=[]

#plot the data
fig2 = plt.figure()
ax2 = fig2.add_subplot(1,1,1)
ax2.scatter(newx, newy, s=10, c='b', marker="s", label='Original Data')
fig2.suptitle("Sunspots")
ax2.set_xlabel("Months from January 1749")
ax2.set_ylabel("Sunspot Number")
plt.legend(loc='upper left');

#deal with the different cases for running average
for x in range(0,newx.size):
	#deals with the case at the beginning of the list
    if(x<5 and newx.size-1>=11):
        newlist=[]
        for l in range(0, 5):
            newlist.append(newy[x])
        newlist.append(newy[x])
        for n in range(x+1,x+6):
            newlist.append(newy[n])
    #deals with the case in the middle of the list
    if(x>=5 and x<newx.size-6):
        newlist=[]
        for l in range(-5, 1):
            newlist.append(newy[x+l])
        for n in range(x+1,x+6):
            newlist.append(newy[n])
    #deals with the case at the end of the list
    if(x>newx.size-6):
        newlist=[]
        for l in range(-5, 1):
            newlist.append(newy[x+l])
        for n in range(x,newx.size-1):
            newlist.append(newy[n+l])
        while(len(newlist)<11):
            newlist.append(newy[x])
    running.append(runaverage(newlist))
#plot the data alongside it's running average
fig2 = plt.figure()
ax1 = fig2.add_subplot(1,1,1)
ax1.scatter(newx, newy, s=10, c='b', marker="s", label='Original Data')
ax1.scatter(newx,running, s=10, c='r', marker="o", label='running average')
fig2.suptitle("Sunspots")
ax1.set_xlabel("Months from January 1749")
ax1.set_ylabel("Sunspot Number")
ax1.legend(loc='upper left')

    
