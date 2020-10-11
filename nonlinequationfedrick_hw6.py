import numpy as np
import matplotlib.pyplot as plt
#method for relaxation method works for both systems of relaxation and single variables x must be a list
def relax(x,f,error):
    newerror=[]
    #create the error vector
    for l in x:
        newerror.append(10)
    i=0
    stop=False
    #interate through function
    while(not stop and i<10000):
        i+=1
        #pass x into a vector function that takes a vector and outputs another vector
        newx=f(x)
        for n in range(0,len(newerror)):
            newerror[n]=abs(x[n]-newx[n])
        x=newx
        c=0
        #checks that all of the values in vector x have the desired error
        for j in newerror:
            if(j<error):
                c+=1
        #if they all have the desired ammount of error we can stop interating
        if(c==len(newerror)):
            stop=True
#this places a limit on the number of times the function can interate
    if(i>10000):
        print("this did not converge")
    return x,i
# newton raphson for a single variable function follows a similar patten as relaxation except for single variables
def NR(f,df,error,x):
    newerror=1
    i=0
    while(error<newerror or i>1000):
        i+=1
        newx=x-(f(x)/df(x))
        print(newx)
        newerror=abs(x-newx)
        x=newx
    if(i>1000):
        print("This did not converge")
    return x
#secant method for a single variable function  follows a similar patten as relaxation except for single variables
def secant(f,error,x1,x2):
    newerror=1
    x=0
    i=0
    while(error<newerror and i<1000):
        i+=1
        try:
            df=(x2-x1)/(f(x2)-f(x1))
        except ZeroDivisionError:
            df=0
        x=x2-f(x2)*df
        x1=x2
        x2=x
        newerror=abs(x1-x2)
    if(i>1000):
        print("This did not converge")
    return x
#this is a help function for the binary search method that produces a plot of the function you are trying to find the root of and allows a manual rechange of the domain
#it will then check if this domain is acceptable
def changerange(a,b,f):
    bad=True
    #print("this is a "+ str(a))
    #print("this is b "+ str(b))
    while(bad):
        y1=f(a)
        y2=f(b)
        #check that the signs of beginning and end of the domain are opposite
        if(((y1<0 and y2)>0) or ((y1>0 and y2)<0)):
            return y1,y2,a,b
        else:
            #if the domain is bad we plot the function with a domain 10X greater and 10X small of the function
            print("this domain is bad please enter a new one")
            x=np.linspace(10*a,10*b,1000)
            y=f(x)
            fig1=plt.figure()
            ax= fig1.add_subplot(2,1,1)
            ax.plot(x,y)
            ax.set_xlabel("X")
            ax.set_ylabel("Y")
            ax.set_title("zoomed out")
            x=np.linspace(0.1*a,0.1*b,1000)
            y=f(x)
            ax2= fig1.add_subplot(2,1,2)
            ax2.plot(x,y)
            ax2.set_xlabel("X")
            ax2.set_ylabel("Y")
            ax2.set_title("zoomed out")
            fig1.tight_layout(pad=3.0)
            plt.show()
            #we request that the user uses the plots produced to change the domain
            a=float(input("Enter the beginning of the domain "))
            b=float(input("Enter the end of the domain "))
#perfom binary search recursivels
def binary(a,b,f,error):
    y1,y2,a,b=changerange(a,b,f)
    newerror=abs(y1-y2)
    avg=(a+b)/2
    yavg=f(avg)
    #if the domain is sufficiently small output a result
    if(newerror<error):
        return avg
    elif(yavg==0):
        return avg
    #splice the based on the value of yavg and the sign of the end points
    elif(yavg<0):
        if(y1<0):
            return binary(avg,b,f,error)
        if(y2<0):
            return binary(a,avg,f,error)
    elif(yavg>0):
        if(y1>0):
            return binary(avg,b,f,error)
        if(y2>0):
            return binary(a,avg,f,error)
#perfomr newton raphson on a system
def NRsys(f,J,error,x):
    newerror=[]
    for l in x:
        newerror.append(10)
    i=0
    stop=False
    while(not stop and i<10000):
        i+=1
        #jacobian matrix become A
        A=J(x)
        # our b matrix is simply our function
        b=f(x)
        #we then solve for change in our function at x
        A=np.array(A)
        b=np.array(b)
        dx=np.linalg.solve(A,b)
        newx=x-dx
        #calculate the new error
        for n in range(0,len(newerror)):
            newerror[n]=abs(x[n]-newx[n])
        c=0
        for j in newerror:
            if(j<error):
                c+=1
        if(c==len(newerror)):
            stop=True
        x=newx
        #print(np.array(x).real)
    if(i>10000):
        print("This did not converge")
    return x,i
