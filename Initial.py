from sympy import *
import numpy as np
import sympy as sp
#from sympy.interactive.printing import init_printing
#init_printing(use_unicode=False, wrap_line=False)
elemNode = 4
nelem = 1 
## Variable
s,t = symbols('s t')

## Initial objects                                       
coords = np.zeros((nelem,1),dtype='object')
D = np.zeros((nelem,1),dtype='object')
E = np.zeros((nelem,1),dtype='object')
nuy = np.zeros((nelem,1),dtype='object')
x = np.zeros((nelem,1),dtype='object')
y = np.zeros((nelem,1),dtype='object')
JacobianMat = np.zeros((nelem,1),dtype='object')
detJ = np.zeros((nelem,1),dtype='object')
B = np.zeros((nelem,1),dtype='object') 
k = np.zeros((nelem,1),dtype='object') 
thickness = np.zeros((nelem,1),dtype='object') 