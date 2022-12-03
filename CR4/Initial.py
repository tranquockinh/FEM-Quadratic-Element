from sympy import *
import numpy as np
import sympy as sp
from sympy.interactive.printing import init_printing
init_printing(use_unicode=False, wrap_line=False)
elemNode = 4
elemEdge = 4
nelem = 2
nnode = 6
nDOF = nnode * 2
K = np.zeros((nDOF,nDOF),dtype='float')
U = np.zeros((nDOF,1),dtype='float')
Fs = np.zeros((nDOF,1),dtype='float')
Fb = np.zeros((nDOF,1),dtype='float')
## Variable
s,t = symbols('s t')
## Initial objects                          
coords = np.zeros((nelem),dtype='object')
D_Store = np.zeros((nelem),dtype='object')
E = np.zeros((nelem),dtype='object')
nuy = np.zeros((nelem),dtype='object')
JacobianMat = np.zeros((nelem),dtype='object')
detJ = np.zeros((nelem),dtype='object')
B_Store = np.zeros((nelem),dtype='object') 
k = np.zeros((nelem),dtype='object') 
u = np.zeros((nelem),dtype='object')
Stress = np.zeros((nelem),dtype='object') 
thickness = np.zeros((nelem),dtype='object')
gamma = np.zeros((nelem),dtype='object')
connect = np.zeros((nelem),dtype='object')
Ts = np.zeros((nelem),dtype='object')
fs = np.zeros((nelem),dtype='object')
fb = np.zeros((nelem),dtype='object')
UnitWeight = np.zeros((nelem),dtype='object')
Ka = np.zeros((nelem),dtype='object')