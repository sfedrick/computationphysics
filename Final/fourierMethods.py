import numpy as np
import matplotlib.pyplot as plt
from pandas.core.common import flatten
#returns The wavelengths from 0 +K1 +K2...KN/2 (-KN/2)-1 -KN -KN-1 -KN-2 ... -K2 -K1
def findK(N):
    K1=[0]
    K2=np.linspace(1,N/2-1,int(N/2)-1).tolist()
    K3=np.linspace(-N/2,-1,int(N/2)).tolist()
    K3[0]=0;
    K=[K1,K2,K3]
    K=list(flatten(K))
    return K
    
def realA(Asmall):
    A1=Asmall[0]
    A2=Asmall[1:-1:1]
    A3=Asmall[:0:-1]
    A=[A1,A2,A3]
    A=list(flatten(A))
    return A
