import numpy as np
import matplotlib.pyplot as plt
from math import exp
from nonlinequation import *
def f(x):
    y=2-exp(-x)-x
    return y
def df(x):
    y=exp(-x)-1
    return y

print(NR(f,df,1e-16,2))
