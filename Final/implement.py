import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
from resonance import *
from fourierMethods import *
def forceinit(Q,size):
    xpoint=0.3*l
    ypoint=0.3*w
    found=False
    for i in range(0,N):
        for j in range(0,N):
            if((XX[i,j]>xpoint-size/2) and (XX[i,j]<xpoint+size/2)):
                if((YY[i,j]>ypoint-size) and (YY[i,j]<ypoint+size)):
                    if(not(found)):
                        Q[i,j]=load
                        found=True
                    else:
                        Q[i,j]=0


N=64
l=1
w=1
E=71*10*9
t=0.003
v=0.3
p=2700
D=(E*t)/(12*(1-v**2))

load=0.5
X=np.linspace(0,w,N)
Y=np.linspace(0,l,N)
YY,XX=np.meshgrid(X,Y)
Q=np.zeros([N,N])
forceinit(Q,l/N)

A=np.zeros([N,N])
x0=np.array([A,A])
a=0
b=1
dx=abs(X[0]-X[1])
dy=abs(Y[0]-Y[1])
dt=0.001
delta=[dt,dx,dy]



plotpointx=0.25
plotpointy=0.87
plotpointx=int(plotpointx/dx)
plotpointy=int(plotpointy/dy)

#X,t=oderes(x0,a,b,N,D,p,Q,delta)
X,t=oderesspec(x0,a,b,N,D,p,Q,delta)
array=list(X.values())
signal=[]
for a in array:
    signal.append(a[plotpointx][plotpointy])


F=np.fft.fft(signal)
K=findK(len(F))

freq = np.fft.fftfreq(len(F), d=dt)

m0=array[2]
mtest=array[5]
m1=array[250]



print(plotpointx)
print(plotpointy)

figure1=plt.figure()
ax=mplot3d.Axes3D(figure1)
ax.plot_surface(XX,YY,m0,cmap='viridis')
ax.set_xlabel("x")
ax.set_ylabel("y")
ax.set_zlabel("Deflection")
figure1.suptitle("Deflection of plate when impulse is first applied")

figure2=plt.figure()
ax=mplot3d.Axes3D(figure2)
ax.plot_surface(XX,YY,m1,cmap='viridis')
ax.set_xlabel("x")
ax.set_ylabel("y")
ax.set_zlabel("Deflection")
figure2.suptitle("Deflection of plate test at 0.25 seconds")


figure3=plt.figure()
ax=figure3.add_subplot()
ax.plot(K,abs(F))
ax.set_xlabel("Wave number K")
ax.set_ylabel("Amplitude")
figure3.suptitle("Time")

figure4=plt.figure()
ax=figure4.add_subplot()
ax.plot(t,signal)
ax.set_xlabel("Time")
ax.set_ylabel("Deflection")
figure4.suptitle("Deflection vs Time")
plt.show()
