import numpy as np
import math
from fourierMethods import *
def derive4(w,dx):
    W=(1/(dx**4))*(w[0]-4*w[1]+6*w[2]-4*w[3]+w[4])
    return W
def mixed4(w,dx,dy):
    #w is an array
    factor=1/((dx**2)*(dy**2))
    base1=4*w[1][1]-2*(w[2,1]+w[1,2]+w[0,1]+w[1,0])
    base2=w[2,2]+w[0,2]+w[0,0]+w[2,0]
    W=factor*(base1+base2)
    return W
def EEhelper1(A,dw,i,j,changerows):
    if(changerows):
        I=range(i-2,i+3)
    else:
        I=range(j-2,j+3)
    N=len(dw)
    e=0
    for n in I:
        if(n<0 or n>=N):
            e=e+1
        elif(changerows):
            dw[e]=A[n,j]
        else:
            dw[e]=A[i,n]

def EEhelper2(A,dw,i,j):
    I=range(i-1,i+2)
    J=range(j-1,j+2)
    Nrow,Ncol=A.shape
    e=-1
    for n in I:
        e+=1
        E=-1
        for s in J:
            if((n<0 or n>=Nrow)):
                break;
            elif(s<0 or s>=Ncol):
                E+=1
            else:
                dw[e,E]=A[n,s]
                E+=1

def EEresonancetrash(input,dx,dy,dt,D,p,q):
    A0=input[0]
    A1=input[1]
    row,col=A1.shape
    A=np.zeros([row,col])
    for i in range(0,row):
        for j in range(0,col):
            dwdx=np.zeros([5])
            dwdy=np.zeros([5])
            dwdxdy=np.zeros([3,3])
            EEhelper1(A1,dwdx,i,j,False)
            EEhelper1(A1,dwdy,i,j,True)
            EEhelper2(A1,dwdxdy,i,j)
            DwDx=derive4(dwdx,dx)
            DwDy=derive4(dwdy,dy)
            DwDxDy=mixed4(dwdxdy,dx,dy)
            timefactor=(((dt**2)*D)/p)
            if(i==0 or i==row-1):
                A[i,j]=0
            elif(j==0 or j==col-1):
                A[i,j]=0
            else:
                A[i,j]=timefactor*(-DwDx-2*DwDxDy-DwDy+(q[i,j]/D))-A0[i,j]+2*A1[i,j]
            if(not math.isfinite(A[i,j])):
                A[i,j]=0
    return A
def checki(A,i,j,di,dj,Ni,Nj):
    i=i+di
    j=j+dj
    if(i<=0 or j<=0):
        return 0
    if(i>=Ni or j>=Nj):
        return 0
    else:
        return A[i,j]

def EEresonance(input,dx,dy,dt,D,p,q):
    A0=input[0]
    A1=input[1]
    row,col=A1.shape

    A=np.zeros([row,col])
    P=(6/(dx**4))+(8/(dx**2)*(dy**2))+(6/(dy**4))
    Q=(-4/(dx**4))-(4/(dx**2)*(dy**2))
    R=(-4/(dy**4))-(4/(dx**2)*(dy**2))
    S=(2/((dx**2)*(dy**2)))
    T=(1/(dx**4))
    U=(1/(dy**4))
    timefactor=-(((dt**2)*D)/p)

    for i in range(0,row):
        for j in range(0,col):

            if(i==0 or i==row-1):
                A[i,j]=0
            elif(j==0 or j==col-1):

                A[i,j]=0
            else:
                #P
                wij=A1[i,j]
                #Q
                wi1j=checki(A1,i,j,1,0,row,col)
                wi_1j=checki(A1,i,j,-1,0,row,col)
                #R
                wij1=checki(A1,i,j,0,1,row,col)
                wij_1=checki(A1,i,j,0,-1,row,col)
                #S
                wi1ji=checki(A1,i,j,1,1,row,col)
                wi_1j1=checki(A1,i,j,-1,1,row,col)
                wi_1j_1=checki(A1,i,j,-1,-1,row,col)
                wi1j_1=checki(A1,i,j,1,-1,row,col)
                #T
                wi2j=checki(A1,i,j,2,0,row,col)
                wi_2j=checki(A1,i,j,-2,0,row,col)
                #U
                wij2=checki(A1,i,j,0,2,row,col)
                wij_2=checki(A1,i,j,0,-2,row,col)
                #dt
                wijk=A1[i,j]
                wijk_1=A0[i,j]
                forceterm=(((dt**2)*q[i,j])/p)
                A[i,j]=timefactor*(P*wij+Q*(wi1j+wi_1j)+R*(wij1+wij_1)\
                +S*(wi1ji+wi_1j1+wi_1j_1+wi1j_1)+T*(wi2j+wi_2j)+U*(wij2+wij_2))\
                +2*wijk-wijk_1+forceterm

    return A



def oderes(x0,a,b,N,D,p,q,delta):
    dt=delta[0]
    dx=delta[1]
    dy=delta[2]
    #x0 is a vector of init values
    #a and b are the beginning and end of the range
    #h is the step size
    #f is a list of function for each ode system
    X=dict()
    j=0
    n=len(x0)
    for i in range(0,n):
        key=str(i)
        X[key]=x0[j]
        j+=1
    t=[a,a]
    end=a
    check1=0
    check2=0
    #how often to check a number from 0 to 1 closer to 1 means less outputs of percentage complete
    PercentCheck=0.1
    created=False
    while(end<=b):
        if(check2-check1>PercentCheck):
            check1=float(end)/b
            print("percent until Euler completes "+str(check1))
        input=[]
        list1=str(len(X)-2)
        list2=str(len(X)-1)
        input.append(X[list1])
        input.append(X[list2])

        #convert the input into np array
        input=np.array(input)
        #sends input to Euler system
        if((end-a)>0.1 and not(created)):
            q=np.zeros(X[list1].shape)
            created=True

        new=EEresonance(input,dx,dy,dt,D,p,q)
        list=str(len(X))
        X[list]=new

        t.append(t[-1]+dt)
        end+=dt
        #check updates the percentage of completion
        check2=float(end)/b
    #return a dictionary of values and a vector of time steps
    return X,t


def Specresonance(input,kx,ky,dt,D,p,q):
    A0=input[0]
    A1=input[1]
    row,col=A1.shape
    timefactor=-(((dt**2)*D)/p)
    A=np.zeros([row,col])
    for i in range(0,row):
        for j in range(0,col):
            (q[i,j]/D) -A0[i,j]+2*A1[i,j]
            w=A1[i,j]
            A[i,j]=timefactor*((-kx[i,j]**4)*w+((-kx[i,j]**2)*(-ky[i,j]**2))*w\
            +(-ky[i,j]**4)*w+ q[i,j]/D)+2*w-A0[i,j]
    return A


def oderesspec(x0,a,b,N,D,p,q,delta):
    dt=delta[0]
    q=np.fft.fft2(q)
    #x0 is a vector of init values
    #a and b are the beginning and end of the range
    #h is the step size
    #f is a list of function for each ode system
    X=dict()
    Xf=dict()
    j=0
    n=len(x0)
    n1=len(x0[j])
    kx=findK(n1)
    ky=findK(n1)
    Ky,Kx=np.meshgrid(kx,ky)
    for i in range(0,n):
        key=str(i)
        X[key]=x0[j]
        Xf[key]=np.fft.fft2(x0[j])
        j+=1
    t=[a,a]
    end=a
    check1=0
    check2=0
    #how often to check a number from 0 to 1 closer to 1 means less outputs of percentage complete
    PercentCheck=0.1
    created=False
    while(end<=b):
        if(check2-check1>PercentCheck):
            check1=float(end)/b
            print("percent until spectral Euler completes "+str(check1))
        input=[]
        list1=str(len(X)-2)
        list2=str(len(X)-1)
        input.append(Xf[list1])
        input.append(Xf[list2])

        #convert the input into np array
        input=np.array(input)
        #sends input to Euler system
        if((end-a)>0.1 and not(created)):
            q=np.zeros(X[list1].shape)
            q=np.fft.fft2(q)
            created=True

        new=Specresonance(input,Kx,Ky,dt,D,p,q)
        list=str(len(X))
        Xf[list]=new
        X[list]=np.fft.ifft2(new)

        t.append(t[-1]+dt)
        end+=dt
        #check updates the percentage of completion
        check2=float(end)/b
    #return a dictionary of values and a vector of time steps
    return X,t
