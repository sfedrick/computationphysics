# -*- coding: utf-8 -*-
"""
Created on Tue Feb 11 23:32:56 2020

@author: shaun
"""
import numpy as np
def firstorder(f1,f2,h):
    answer=(f1-f2)/h
    return answer
def firstorderwhole(listx,listy,h):
    L=len(listy)
    y=[]
    x=[]
    for i in range(0,L-1):
        y.append(firstorder(listy[i],listy[i+1],h))
        x.append(listx[i])
    return [x,y]

def secondorder(f1,f3,h):
    y=[]
    x=[]
    answer=(f3-f1)/(2*h)
    return answer
def secondorderwhole(listx,listy):
    L=len(listy)
    y=[]
    x=[]
    for i in range(1,L-1):
        h=abs(listx[i]-listx[i-1])
        y.append(secondorder(listy[i-1],listy[i],h))
        x.append(listx[i])
    x=np.array(x)
    y=np.array(y)
    return x,y
def fourthorder(f1,f2,f3,f4,h):
    f=(f1-8*f2+8*f3-f4)/12*h
    return f
def fourthorderwhole(listx,listy,h):
    y=[]
    x=[]
    L=len(listy)
    for i in range(3,L-2,3):
        y.append(fourthorder(listy[i-2],listy[i-1],listy[i+1],listy[i+2],h))
        x.append(listx[i])
    return x,y

def fourthorders(f1,f2,f3,f4,h):
    A=np.matrix([[1,1,1,1],[0,1,2,3],[0,0.5,2,9/2],[0,1/6,8/6,27/6]],float)
    B=np.matrix([[0],[-1/h],[0],[0]],float)
    A_inverse= np.linalg.inv(A)
    a=A_inverse*B
    answer=a[0]*f1+a[1]*f2+a[2]*f3+a[3]*f4
    return answer

def tridiag(a, b, c, k1=-1, k2=0, k3=1):
    return np.diag(a, k1) + np.diag(b, k2) + np.diag(c, k3)


#pade scheme
def pade(a,b,c,d,l,F):
    A=[]
    B=[]
    C=[]
    D=[]
    B.append(b[0])
    C.append(c[0])
    print("this is l length:"+ str(l))
    for k in range(0,l-2):
        A.append(a[1])
        B.append(b[1])
        C.append(c[1])
    B.append(b[2])
    A.append(a[2])
    D.append(d[0][0]*F[0]+d[0][1]*F[1]+d[0][2]*F[2])
    print("this is F length:"+str(len(F)))
    for k in range(0,len(F)-2):
        D.append(d[1][0]*(F[k+2]-F[k]))
    D.append(d[2][0]*F[-3]+d[2][1]*F[-2]+d[2][2]*F[-1])

    T = tridiag(A, B, C)
    print(T)
    alpha=[B[0]]
    beta=[]
    victor=[D[0]]
    N=len(A)
    f=[]
    for k in range(0,N):
        beta.append((A[k]/alpha[k]))
        alpha.append(B[k+1]-C[k]*beta[k])
    for k in range(0,N):
        victor.append(D[k+1]-beta[k]*victor[-1])
    f.append(victor[-1]/alpha[-1])
    for k in range(N-1,-1,-1):
       f.append((victor[k]-C[k]*f[-1])/alpha[k])
    f.reverse()
    return f

#use this if you do not have a triangular structure
def badpade(a,b,c,d,l,F):
    A=[]
    B=[]
    C=[]
    D=[]
    B.append(b[0])
    C.append(c[0])
    for k in range(0,l-2):
        A.append(a[1])
        B.append(b[1])
        C.append(c[1])
    B.append(b[2])
    A.append(a[2])

    D.append([d[1][0]*(F[1]-F[l-1])])
    for k in range(0,len(F)-2):
        D.append([d[1][0]*(F[k+2]-F[k])])
    D.append([d[1][0]*(F[0]-F[l-1])])
    T = tridiag(A, B, C)
    T[0,2]=c[2]
    T[-1,-3]=a[0]
    D=np.array(D)
    IT=np.linalg.inv(T)*(D.transpose())
    f=np.matmul(IT,D)
    return f
d=np.array([[1,2,3],[3,0,0],[7,8,9]], int)
ans=badpade([-5,1,2],[1,4,1],[2,1,82],d,5,[1,1,1,1,1])
#print(pade([0,1,2],[1,4,1],[2,1,0],d,5,[1,1,1,1,1]))
print(ans)
