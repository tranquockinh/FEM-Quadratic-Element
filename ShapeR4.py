import numpy as np
import sympy as sp
from sympy import *

# Predefine
nelem = 1
nodeE = 4
x,y = symbols('x y')
r,s = symbols('r s')
coords = np.zeros((nelem,1),dtype = 'object')
N = np.zeros((nelem,1),dtype = 'object')
xmap = np.zeros((nelem,1),dtype = 'object')
ymap = np.zeros((nelem,1),dtype = 'object')
J = np.zeros((nelem,1),dtype = 'object')
detJ = np.zeros((nelem,1),dtype = 'object')
G = np.zeros((nelem,1),dtype = 'object')
D = np.zeros((nelem,1),dtype = 'object')
E = np.zeros((nelem,1),dtype = 'object')
nuy = np.zeros((nelem,1),dtype = 'object')
t = np.zeros((nelem,1),dtype = 'object')
B = np.zeros((nelem,1),dtype = 'object')
SF = np.zeros((nelem,1),dtype = 'object')
Binterpo = np.zeros((nelem,1),dtype = 'object')
K = np.zeros((nelem,1),dtype = 'object')

# Shape function

InterpoMatrix = 1/4*sp.Matrix([[1,-s,-r,r*s],
                               [1,-s,r,-r*s],
                               [1,s,r,r*s],
                               [1,s,-r,-r*s]])
               
# Coordinate  
E[0,0] = 30E06
nuy[0,0] = 0.3   
t[0,0] = 1         
coords[0,0] = np.array([[1,0],
                        [2,0],
                        [2.25,1.5],
                        [1.25,1]],dtype='float')
rscoords = np.array([[-1,-1],
                     [+1,-1],
                     [+1,+1],
                     [-1,+1]],dtype='float')   
a = 1/2*np.abs((rscoords[0,0] - rscoords[1,0]))
b = 1/2*np.abs((rscoords[0,1] - rscoords[-1,1]))

# Gausian interpolation points  
r_GausianPoint = np.array([1/np.sqrt(3),-1/np.sqrt(3)])
s_GausianPoint = np.array([1/np.sqrt(3),-1/np.sqrt(3)])
W1 = 1
W2 = 1
for i in range(nelem):
    N[i,0] = np.zeros((1,nodeE),dtype='object')
    J[i,0] = np.zeros((2,2),dtype='object')
    SF[i,0] = np.zeros((2,8),dtype='object')
    B[i,0] = np.zeros((3,8),dtype='object')
    
    '''
    N[i,0][0,0] = 1/4 * (1-r)*(1-s)
    N[i,0][0,1] = 1/4 * (1+r)*(1-s)
    N[i,0][0,2] = 1/4 * (1+r)*(1+s)
    N[i,0][0,3] = 1/4 * (1-r)*(1+s)
    '''
    N[i,0][0,0] = np.sum(InterpoMatrix[0,:])
    N[i,0][0,1] = np.sum(InterpoMatrix[1,:])
    N[i,0][0,2] = np.sum(InterpoMatrix[2,:])
    N[i,0][0,3] = np.sum(InterpoMatrix[3,:])
    
    N1,N2,N3,N4 = N[i,0][0,0],N[i,0][0,1],N[i,0][0,2],N[i,0][0,3]
    
    SF[i,0][0,:4],SF[i,0][1,4:]  = N[i,0], N[i,0]
    
    for j in range(8):
        B[i,0][0,j] = 1/a*diff(SF[i,0][0,j],r)
        B[i,0][1,j] = 1/b*diff(SF[i,0][1,j],s)
        B[i,0][2,:4],B[i,0][2,4:] = B[i,0][1,4:],B[i,0][0,:4]
    
    xmap[i,0] = np.matmul(N[i,0], coords[i,0][:,0])
    ymap[i,0] = np.matmul(N[i,0], coords[i,0][:,1])
    
    J[i,0][0,0] = diff(xmap[i,0][0],r)
    J[i,0][0,1] = diff(ymap[i,0][0],r)
    J[i,0][1,0] = diff(xmap[i,0][0],s)
    J[i,0][1,1] = diff(ymap[i,0][0],s)
    
    J11,J12,J21,J22 = J[i,0][0,0],J[i,0][0,1],J[i,0][1,0],J[i,0][1,1]
    
    ##J11,J12,J21,J22 = 1/2,1/4*(0.5-0.5*s),0.5,1/4*(2.5-0.5*r)
    ## Convert to symbolic expression
    ##J[i,0] = sp.Matrix(J[i,0])
    ##detJ[i,0] = (J[i,0]).det()
    
    detJ[i,0] = J11*J22 - J12*J21
    print(detJ)
    
    ## Geometric matrix
    G[i,0] = 1/detJ[i,0] * sp.Matrix([[J22,-J22,0,0],
                                      [0,0,-J21,J11],
                                      [-J21,J11,J22,-J12]])
                                 
    D[i,0] = (E[i,0]/(1-nuy[i,0]**2)) * np.array([[1, nuy[i,0], 0],
                                                     [nuy[i,0], 1, 0],
                                                     [0, 0, (1-nuy[i,0])/2]])
    '''
    D[i,0] = (E[i,0]/((1+nuy[i,0])*(1-2*nuy[i,0]))) * np.array([[1-nuy[i,0], nuy[i,0], 0],
                                                                   [nuy[i,0], 1-nuy[i,0], 0],
                                                                   [0, 0, (1-2*nuy[i,0])/2]])    
    '''                                                               
for i in range(nelem):
    Binterpo[i,0] = np.zeros((3,8),dtype='object')
    K[i,0] = np.zeros((nodeE*2,nodeE*2),dtype='object')
    
    for k in range(2):
        for l in range(2):
            Binterpo[i,0] = sp.Matrix(Binterpo[i,0])
            B[i,0] = sp.Matrix(B[i,0])
            Binterpo[i,0] = B[i,0].subs([(r,r_GausianPoint[k]),(s,s_GausianPoint[l])])
            detJ[i,0] = detJ[i,0].subs([(r,r_GausianPoint[k]),(s,s_GausianPoint[l])])
            K[i,0] += W1*W2*t[i,0]*np.matmul(np.matmul(Binterpo[i,0].T,D[i,0]),Binterpo[i,0])*detJ[i,0]
            print(r_GausianPoint[k],s_GausianPoint[l])
            
print(D[0,0])
print()
print(J)
print()
print(B[0,0])
print()
print(K[0,0])

