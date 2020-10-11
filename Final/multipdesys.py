import numpy as np
#perform an expilicit march using euler approximation
def eulerstep(input,grid,t,f,delta):
    new=[]
    #get shape of input
    row,col=input[0].shape
    #create a new array to contain the new time step
    for i in range(0,len(input)):
        array=np.zeros([row,col],float)
        new.append(array)
    #for each element in the newarrays we need a value at x,y
    #or i j in matrix notation
    for i in range(0,row):
        for j in range(0,col):
            n=0
            X=[i,j]
            for matrix in new:
                matrix[i,j]=f[n](X,grid,input,delta)
                n+=1
    return new

def eulerstep1d(input,grid,t,f,delta):
    new=[]
    #get shape of input
    row=len(input[0])
    #create a new array to contain the new time step
    for i in range(0,len(input)):
        array=np.zeros([row],float)
        new.append(array)
    #for each element in the newarrays we need a value at x,y
    #or i j in matrix notation
    for i in range(0,row):
            n=0
            X=[i]
            for matrix in new:
                matrix[i]=f[n](X,grid,input,delta)
                n+=1
    return new
def multipdesys1D(Lm,grid,f,a,b,dt):
    t=[]
    #find the difference in the grid in the x and y direction
    dx=abs(grid[0][0]-grid[0][1])
    #creates a list of deltas so we can easily call dt,dx, and dy
    delta=[dt,dx]
    #create a dictionary that will contain
    # a list of arrays for each variable in the odesolversystem
    #for example if we had a coupled system u and v
    # X contain two list one for u and one for v
    # each element in these lists is a matrix containing all the x,y positions
    #where the last element in the list corresponds to the final time step
    X=dict()
    #create a index for number of coupled systems
    n=len(Lm)
    for i in range(0,n):
        #create a key for each system in the coupled system
        key=str(i)
        X[key]=[Lm[i]]
    check=0
    while a<b:
        #add to t list which is the list of time steps
        t.append(a)
        input=[]
        new=[]
        #pulls the last time step from each list in the dictionary
        for list in X:
            #this creates the array that is fed into the Euler system algorithm
            input.append(X[list][-1])
        #performs euler step and returns a matrix for each
        #variable in the coupled system
        new=eulerstep1d(input,grid,t[-1],f,delta)
        i=0
        #for every list in x we append the new time marhc
        for list in X:
            #appends the output to the list
            X[list].append(new[i])
            i=i+1
        a=a+dt
        check+=1
    # X will be a dictionary of lists each list will contain arrays corresponding
    # to the 2d grid that progressed in time
    return X,t

#solve a pde system in 2 dimensions
def multipdesys2dEE(Lm,grid,f,a,b,dt):
    t=[]
    #find the difference in the grid in the x and y direction
    dx=abs(grid[0][0,0]-grid[0][1,0])
    dy=abs(grid[1][0,0]-grid[1][0,1])
    #creates a list of deltas so we can easily call dt,dx, and dy
    delta=[dt,dx,dy]
    #create a dictionary that will contain
    # a list of arrays for each variable in the odesolversystem
    #for example if we had a coupled system u and v
    # X contain two list one for u and one for v
    # each element in these lists is a matrix containing all the x,y positions
    #where the last element in the list corresponds to the final time step
    X=dict()
    #create a index for number of coupled systems
    n=len(Lm)
    for i in range(0,n):
        #create a key for each system in the coupled system
        key=str(i)
        X[key]=[Lm[i]]
    check=0
    while a<b:
        #add to t list which is the list of time steps
        t.append(a)
        input=[]
        new=[]
        #pulls the last time step from each list in the dictionary
        for list in X:
            #this creates the array that is fed into the Euler system algorithm
            input.append(X[list][-1])
        #performs euler step and returns a matrix for each
        #variable in the coupled system
        new=eulerstep(input,grid,t[-1],f,delta)
        i=0
        #for every list in x we append the new time marhc
        for list in X:
            #appends the output to the list
            X[list].append(new[i])
            i=i+1
        a=a+dt
        check+=1
    # X will be a dictionary of lists each list will contain arrays corresponding
    # to the 2d grid that progressed in time
    return X,t

#a is the lower diagnol
#b is the diagnol
#c is the upper diagnol
#d are the constants being dolved for
def TDMAsolver(a, b, c, d):
    '''
    TDMA solver, a b c d can be NumPy array type or Python list type.
    refer to http://en.wikipedia.org/wiki/Tridiagonal_matrix_algorithm
    and to http://www.cfd-online.com/Wiki/Tridiagonal_matrix_algorithm_-_TDMA_(Thomas_algorithm)
    '''
    nf = len(d) # number of equations
    ac, bc, cc, dc = map(np.array, (a, b, c, d)) # copy arrays
    for it in range(1, nf):
        mc = ac[it-1]/bc[it-1]
        bc[it] = bc[it] - mc*cc[it-1]
        dc[it] = dc[it] - mc*dc[it-1]

    xc = bc
    xc[-1] = dc[-1]/bc[-1]

    for il in range(nf-2, -1, -1):
        xc[il] = (dc[il]-cc[il]*xc[il+1])/bc[il]

    return xc
def implicithelp(input,grid,t,f,fI,delta):
    row,col=input[0].shape
    dmatrix=[]
    for i in range(0,len(input)):
        array=np.zeros([row,col],float)
        dmatrix.append(array)
    for i in range(0,row):
        for j in range(0,col):
            n=0
            X=[i,j]
            for matrix in dmatrix:
                matrix[i,j]=f[n](X,grid,input,delta)
                n+=1
    return dmatrix
def diag(input,grid,delta,f):
    row,col=input[0].shape
    A=[]
    B=[]
    C=[]
    i=0
    for matrix in input:
        a,b,c=f[i](row,grid,delta)
        a=np.array(a)
        b=np.array(b)
        c=np.array(c)
        A.append(a)
        B.append(b)
        C.append(c)
        i+=1
    return [A,B,C]
def impliciteulerstepT(input,grid,t,f,g,fIx,fIy,delta,flip):

    new=[]
    #create a tridiagnol matrix for for x march and y march
    #a is the lower diagnol
    #b is the diagnol
    #c is the upper diagnol
    ax,bx,cx=diag(input,grid,[delta[0],delta[1]],fIx)
    ay,by,cy=diag(input,grid,[delta[0],delta[2]],fIy)
    row,col=input[0].shape
    #we create an array of zeros to fill with new information
    for i in range(0,len(input)):
        array=np.zeros([row,col],float)
        new.append(array)
    #check in which order to perform the explicit and implicit march
    if(flip%2!=0):

        dmatrix1=implicithelp(input,grid,t,f,fIy,delta)
        for j in range(0,col):
            n=0
            for matrix in new:
                ans=TDMAsolver(ax[n],bx[n],cx[n],dmatrix1[n][:,j])
                matrix[:,j]=ans
                n+=1
        dmatrix2=implicithelp(new,grid,t,g,fIx,delta)
        for i in range(0,row):
            n=0
            for matrix in new:
                ans=TDMAsolver(ay[n],by[n],cy[n],dmatrix2[n][i,:])
                matrix[i,:]=ans
                n+=1
        return new
    #check in which order to perform the explicit and implicit march
    else:
        dmatrix1=implicithelp(input,grid,t,g,fIx,delta)
        for i in range(0,row):
            n=0
            for matrix in new:
                ans=TDMAsolver(ay[n],by[n],cy[n],dmatrix1[n][i,:])
                matrix[i,:]=ans
                n+=1
        dmatrix2=implicithelp(new,grid,t,f,fIy,delta)
        for j in range(0,col):
            n=0
            for matrix in new:
                ans=TDMAsolver(ax[n],bx[n],cx[n],dmatrix2[n][:,j])
                matrix[:,j]=ans
                n+=1
        return new
#f is the right hand side when going in the x direction in adi
#g is the right hand side when going in the y direction in adi
#fIx is the diagnol or left hand side when going in the x direction
#fIy is the diagnol or left hand side going in the y direction
def multipdesys2dIET(Lm,grid,f,g,fIx,fIy,a,b,dt):
    t=[]
    flip=0
    #find the difference in the grid in the x and y direction
    dx=abs(grid[0][0,0]-grid[0][1,0])
    dy=abs(grid[1][0,0]-grid[1][0,1])
    delta=[dt/2,dx,dy]
    X=dict()
    n=len(Lm)
    for i in range(0,n):
        key=str(i)
        X[key]=[Lm[i]]
    check=0
    while a<b:
        flip+=1
        t.append(a)
        input=[]
        new=[]
        for list in X:
            #this creates the array that is fed into the Euler system algorithm
            input.append(X[list][-1])
        new=impliciteulerstepT(input,grid,t[-1],f,g,fIx,fIy,delta,flip)
        if(check==0):
            print("this is n ="+ str(check))
            print(input)
        if(check>10 and check<20 ):
            print("this is n ="+ str(check))
            print(new)
        i=0
        for list in X:
            #appends the output to the list
            X[list].append(new[i])
            i=i+1
        a=a+dt
        check+=1
    return X,t
