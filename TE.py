## One element
import numpy as np
from sympy import *
from sympy.interactive.printing import init_printing
init_printing(use_unicode=False, wrap_line=False)

# Geometry and initialization
x, y = symbols('x y')
coords = Matrix([[0, 0, -1],
                 [1, 2, 0],
                 [2, 0, 1]])
## Natural coordinate matrix
theta, theta[:,0], theta[:,1:] = np.zeros((3,3)), np.array([1,1,1]), coords[:,1:]
A = 1/2 * np.linalg.det(theta)

# Materials
t = 1 ## thickness
Ed = 30E6 ## Young modulus
nuy = 0.25 ## Poisson's ratio

# If soil material and soil structure interation
##phi = 38 * np.pi/180 ## if Soil or sand --- internal friction angle
# Material physical properties
##gammaDam = 22
##gammaSand = 16
##Ka = np.tan(np.pi/4 - phi/2)**22 ## Earth pressure at static
##SurLoad = 50
## Computed from above input
##Gd = Ed * (2*(1+nuy))
##Lambda = (nuy*Ed)/((1-2*nuy) * (1+nuy))

D = (Ed/(1-nuy**2)) * np.array([[1, nuy, 0],
                               [nuy, 1, 0],
                               [0, 0, (1-nuy)/2]])
# Shape matrix
natural_disp = eye(3)
#a = np.zeros((3,3))
#for i in range(3):
a = np.matmul(np.linalg.inv(theta), natural_disp)
u = Matrix([1, x, y]).T * a
N = np.zeros((2,6),dtype=object)
N[0,:3], N[1,3:] = u, u
B = np.zeros((2,6))
for i in range(len(N[0,:])):
    B[0,i] = diff(N[0,i],x) 
    B[1,i] = diff(N[1,i],y)
B = np.reshape(B.T, (2,6))
B = np.concatenate((B, np.zeros((1,6))), axis=0)
B[2,0::2], B[2,1::2] = B[1,1::2], B[0,0::2]
print()

# Stiffness matrix
K = (t * A) * np.matmul(np.matmul(B.T, D), B)
print(B)

# Compute stress if displacement known 
u = np.array([0, 0.0025, 0.0012, 0, 0, 0.0025])[:, np.newaxis]
print(u)
S = np.matmul(np.matmul(D, B), u)
print('Stresses:Sx, Sy and Txy\n', S)

# Pricipal stress and principal angle
S1 = (S[0,0]+S[1,0])/2 + np.sqrt(((S[0,0]-S[1,0])/2)**2 + (S[2,0])**2)
S2 = (S[0,0]+S[1,0])/2 - np.sqrt(((S[0,0]-S[1,0])/2)**2 + (S[2,0])**2)
pricipal_ang = np.degrees(1/2 * np.arctan(np.divide(2*S[2,0], (S[0,0]-S[1,0]))))
print()
print('S1:\n{}\nS2:\n{}\nPrincipal angle:\n{}'.format(S1, S2, pricipal_ang))