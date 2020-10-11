import pandas as pd
import  scipy as sp
import matplotlib.pyplot as plt
from io import StringIO
import numpy as np
from scipy.interpolate import interp1d
from derivativesfedrick_hw6 import *
from integrationfedrick_hw6 import *
from nonlinequationfedrick_hw6 import *
#input data
data=np.loadtxt('cambfedrick_hw6.dat')
df = pd.DataFrame(data)
df.to_excel("tester.xls", index=False)
l=3

n=len(data[:,0])
#pulls 200 points from the cmb
x200=data[0:n:10,0]
y200=data[0:n:10,1]
#perform an interpolation on the 200 data points
cf=interp1d(x200,y200,kind='cubic')
xnew=np.linspace(2,2191)
#grabs the whole data set to compare the error
xwhole=data[:,0]
ywhole=data[:,1]
yerror=abs(cf(xnew)-ywhole[0:len(xnew)])/ywhole[0:len(xnew)]
yerror=np.average(yerror)
print("this is the total error "+ str(yerror))

#compute more points between the orginal 2000 points
sx=np.linspace(2,2191,10000)
sy=cf(sx)
#convert into the desired integral form
cy=((cf(sx))*2*(np.pi))/(l*(l+1))
h=sx[0]-sx[1]
#perfomr fourth order derivative
dfxc,dfyc=fourthorderwhole(sx,cy,h)
#compute intergral over all the points
#define integrand
def cI(x):
    y=(((cf(x)*2*(np.pi))/(l*(l+1)))*(l/(2*np.pi)))
    return y
Fxc,Fyc,sum=twhole(10000,2,2191,cI)
fig1=plt.figure()

#make plots
ax=fig1.add_subplot()
ax.scatter(x200,y200, label="200 points from cmb")
ax.plot(sx,sy, label="Cubic fit of 200 points")
ax.set_title(" Power law of CMB ")
ax.set_xlabel("Power")
ax.set_ylabel("C_{TT}")
ax.legend()

fig2=plt.figure()
ax2=fig2.add_subplot()
ax2.set_title("Derivative of CMB")
ax2.plot(dfxc,dfyc)

fig3=plt.figure()
ax3=fig3.add_subplot()
ax3.plot(Fxc,Fyc)
ax3.set_title("Integral of CMB")
plt.show()
