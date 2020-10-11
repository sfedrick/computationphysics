import numpy as np
#create matrix of voltage equations from the circuit
a=np.array([
[4,-1,-1,-1],
[-1,3,0,-1],
[-1,0,3,-1],
[-1,-1,-1,4]])
#set constants
b=np.array([5,0,5,0])\
#solve matrix
x=np.linalg.solve(a,b)
print(x)
