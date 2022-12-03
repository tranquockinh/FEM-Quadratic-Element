from sympy import *
import numpy as np
from sympy.interactive.printing import init_printing
init_printing(use_unicode=False, wrap_line=False)
elemNode = 4
x,y = symbols('x y')
nelem = 1                                            ## Numnber of elements
nnode = nelem * elemNode                                       ## Numnber of nodes
nDOF = 2 * nnode                                        ## Numnber DOFs
coords = np.empty((nelem,1), dtype='object')            ## Element coordinates
theta = np.empty((nelem,1), dtype='object')             ## Natural coordinate
A = np.empty((nelem,1), dtype='object')                 ## Element area
t = np.empty((nelem,1), dtype='object')                 ## Thickness
E = np.empty((nelem,1), dtype='object')                 ## Elastic modulus
nuy = np.empty((nelem,1), dtype='object')               ## Possion's rato
D = np.empty((nelem,1), dtype='object')                 ## Material elaticity matrix
abcd = np.empty((nelem,1), dtype='object')               ## alpha, beta, gamma
Lu = np.empty((nelem,1), dtype='object')                ## Displacement operator
u = np.empty((nelem,1), dtype='object')                 ## Element displacement  
N = np.empty((nelem,1), dtype='object')                 ## Linear operator
B = np.empty((nelem,1), dtype='object')                 ## Shape function
k = np.empty((nelem,1), dtype='object')                 ## Stiffness matrix
K = np.empty((nelem,1), dtype='object')                 ## Global stiffness matrix
connect = np.empty((nelem,1), dtype='object')           ## Element connectivity
DOF = np.empty((nelem,1), dtype='object')               ## Global degree of freedom
gamma = np.empty((nelem,1), dtype='object')             ## Element material density
c = np.empty((nelem,1), dtype='object')                 ## Element lateral gravitation coefficient
bdloads = np.empty((nelem,1), dtype='object')           ## Body loads       
cloads = np.empty((nelem,1), dtype='object')            ## Concentrated loads
sfloads = np.empty((nelem,1), dtype='object')           ## equivalent concentrated loads
F = np.empty((nelem,1), dtype='object')                 ## Global force vector
Stress = np.empty((nelem,1), dtype='object')            ## Element stresses 
Princ_Stress =  np.empty((nelem,1), dtype='object')
kint =  np.empty((nelem,1), dtype='object')
coords[0,0] = np.array([[0,0],
                        [8,0],
                        [8,4],
                        [0,4]], dtype='float') ## Element 1
connect[0,0] = np.array([0,1,2,3])
t[0,0] = 1
E[0,0] = 30E6
nuy[0,0] = 0.3
  
def local_stiffness():
    for i in range(nelem):
        kint[i,0] = np.zeros((8,8),dtype='object') 
        theta[i,0], theta[i,0][:,0] = np.zeros((4,4)), [1,1,1,1]
        theta[i,0][:,1:3] = coords[i,0][:,0:]
        theta[i,0][:,-1] = coords[i,0][:,0] * coords[i,0][:,1]
        
        D[i,0] = (E[i,0]/(1-nuy[i,0]**2)) * np.array([[1, nuy[i,0], 0],
                                                     [nuy[i,0], 1, 0],
                                                     [0, 0, (1-nuy[i,0])/2]])
        '''
        D[i,0] = (E[i,0]/((1+nuy[i,0])*(1-2*nuy[i,0]))) * np.array([[1-nuy[i,0], nuy[i,0], 0],
                                                                   [nuy[i,0], 1-nuy[i,0], 0],
                                                                   [0, 0, (1-2*nuy[i,0])/2]])
        '''
        abcd[i,0] = np.matmul(np.linalg.inv(theta[i,0]), eye(4))
        ## Compute general displacement operator
        Lu[i,0] = np.matmul(Matrix([1,x,y,x*y]).T, abcd[i,0]) 
        N[i,0],N[i,0][0,0::2],N[i,0][1,1::2] = np.zeros((2,8),dtype='object'),Lu[i,0],Lu[i,0]
        ## Computation of shape function
        B[i,0] = np.zeros((3,8),dtype='object')
        for j in range(8):
            B[i,0][0,j],B[i,0][1,j] = diff(N[i,0][0,j],x), diff(N[i,0][1,j],y)
            B[i,0][2,0::2], B[i,0][2,1::2] = B[i,0][1,1::2], B[i,0][0,0::2] 
        ## Construction of element stiffness matrix    
        k[i,0] = np.matmul(np.matmul(B[i,0].T, D[i,0]), B[i,0])
         
        for m in range(8):
            for l in range(8):
                l1 = coords[i,0][0,0] 
                l2 = coords[i,0][1,0] - l1
                l3 = coords[i,0][1,1]
                l4 = coords[i,0][2,1] - l3
                k[i,0][m,l] = round(integrate(integrate(k[i,0][m,l]*t[i,0],(x,l1,l2)),(y,l3,l4)),3)
                
    return k,N,B,D
    
k,N,B,D = local_stiffness()
print(B[0,0])
print(k[0,0])