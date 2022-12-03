import sympy as sp
import numpy as np
from sympy import *
t,s = sp.symbols('s t')

func = np.array([1,s,t,s*t,s**2,t**2,s**2*t,s*t**2])
MatrixFunc = np.zeros((8,8),dtype='object')

for i in range(8):
    MatrixFunc[i,:] = func

b = np.eye(8,8)
MatrixFunc = sp.Matrix(MatrixFunc)
MatrixFuncinv = MatrixFunc.inv()