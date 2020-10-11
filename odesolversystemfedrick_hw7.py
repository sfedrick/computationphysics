import numpy as np
from scipy.optimize import brentq
# finds the k values for rk4
def rk4step(yn,tn,f,h):
    #note that k1 k2 k3 and k4 are lists the length of these lists
    #corresponds to the number of odes in your system
    a=1/2.0
    k1=h*rk4help(yn,tn,f)
    k2=h*rk4help(yn+a*k1,tn+h/2.0,f)
    k3=h*rk4help(yn+a*k2,tn+h/2.0,f)
    k4=h*rk4help(yn+k3,tn+h,f)
    yn1=yn+(1/6.0)*k1+(1/3.0)*(k2+k3)+(1/6.0)*k4
    return yn1
def rk4help(yn,tn,f):
    #helper function that finds the k value from each function and appends it to K
    k=[]
    for func in f:
        k.append(func(yn,tn))
    k=np.array(k)
    return k
# take in a system of odes
def odesysRk4(x0,a,b,h,f):
    #x0 is a vector of init values
    #a and b are the beginning and end of the range
    #h is the step size
    #f is a list of function for each ode system
    X=dict()
    j=0
    n=len(x0)
    for i in range(0,n):
        key=str(i)
        X[key]=[x0[j]]
        j+=1
    t=[a]
    end=a
    check1=0
    check2=0
    #how often to check a number from 0 to 1 closer to 1 means less outputs of percentage complete
    PercentCheck=0.01
    while(end<b):
        if(check2-check1>PercentCheck):
            check1=float(end)/b
            print("percent until Rk4 completes "+str(check1))
        input=[]
        for list in X:
            #this creates the array that is fed into the rk4 algorithm
            input.append(X[list][-1])
        #convert the input into np array
        input=np.array(input)
        #sends input to rk4 algorithm
        new=rk4step(input,t[-1],f,h)
        i=0
        for list in X:
            #appends the output to the list
            X[list].append(new[i])
            i=i+1
        t.append(t[-1]+h)
        end+=h
        #check updates the percentage of completion
        check2=float(end)/b
    return X,t
#uses shooting method to solve boundary value problems
def linshot(g1,g2,L,T,h,f):
    T1,t1=odesysRk4(g1,0,L,h,f,2)
    T2,t1=odesysRk4(g2,0,L,h,f,2)
    G1,dydx1=T1.values()
    G2,dydx2=T2.values()
    row1=[1,1]
    row2=[G1[-1],G2[-1]]
    A=np.array([row1,row2])
    IA=np.linalg.inv(A)
    C=np.dot(IA,[1,T])
    s=np.multiply(C[0],G1)+np.multiply(C[1],G2)
    return s,t1

def tridiag(a, b, c, k1=-1, k2=0, k3=1):
    return np.diag(a, k1) + np.diag(b, k2) + np.diag(c, k3)


#pade scheme
def direct(a,b,c,d,l,R):
    A=[]
    B=[]
    C=[]
    D=[]
    t=np.linspace(R[0],R[1],l)
    h=abs(t[0]-t[1])
    B.append(b(t[0],h,1))
    C.append(c(t[0],h,1))
    for k in range(1,l-1):
        A.append(a(t[k],h,0))
        B.append(b(t[k],h,0))
        C.append(c(t[k],h,0))
    B.append(b(t[-1],h,-1))
    A.append(a(t[-1],h,-1))
    D.append(d(t[0],h,1))
    for k in range(1,l-1):
        D.append(d(t[k],h,0))
    D.append(d(t[k],h,-1))

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
    return f,t
def eIlinstep(yn,tn,f,df,h):
    l=len(f)
    output=[]
    x=[]
    I=np.array([[1,0,0],[0,1,0],[0,0,1]])
    B=I-(h/2)*df(yn)

    IB=np.linalg.inv(B)

    for i in range(0,l):
        x.append(h*f[i](yn,tn+h))
    yn=np.transpose(yn)
    output=yn+h*np.dot(IB,x)
    output=np.transpose(output)
    return output

def eulerIlin(x0,a,b,h,f,df,n):
    X=dict()
    j=0
    for i in range(0,n):
        key=str(i)
        X[key]=[x0[j]]
        j+=1
    t=[a]
    end=a
    check1=0
    check2=0
    while(end<b):
        if(check2-check1>0.01):
            check1=float(end)/b
            print("percent until Linear Trapezoidal completes "+str(check1))
        input=[]
        for list in X:
            input.append(X[list][-1])
        input=np.array(input)
        newy=eIlinstep(input,t[-1],f,df,h)
        i=0
        for list in X:
            X[list].append(newy[i])
            i=i+1
        t.append(t[-1]+h)
        end+=h
        check2=float(end)/b
    return X,t
