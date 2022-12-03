from sympy import *
import numpy as np
from sympy.interactive.printing import init_printing
init_printing(use_unicode=False, wrap_line=False)

# Coordinate definition and initializations
x,y = symbols('x y')
nelem = 2                                               ## Numnber of elements
nnode = 2 * nelem                                       ## Numnber of nodes
coords = np.empty((nelem,1), dtype='object')            ## Element coordinates
theta = np.empty((nelem,1), dtype='object')             ## Natural coordinate
A = np.empty((nelem,1), dtype='object')                 ## Element area
t = np.empty((nelem,1), dtype='object')                 ## Thickness
E = np.empty((nelem,1), dtype='object')                 ## Elastic modulus
nuy = np.empty((nelem,1), dtype='object')               ## Possion's rato
D = np.empty((nelem,1), dtype='object')                 ## Material elaticity matrix
abg = np.empty((nelem,1), dtype='object')               ## alpha, beta, gamma
u = np.empty((nelem,1), dtype='object')                 ## Displacement
N = np.empty((nelem,1), dtype='object')                 ## Linear operator
B = np.empty((nelem,1), dtype='object')                 ## Shape function
k = np.empty((nelem,1), dtype='object')                 ## Stiffness matrix
K = np.empty((nelem,1), dtype='object')                 ## Global stiffness matrix
connect = np.empty((nelem,1), dtype='object')           ## Element connectivity
DOF = np.empty((nelem,1), dtype='object')               ## Global degree of freedom
gamma = np.empty((nelem,1), dtype='object')             ## Element material density
c = np.empty((nelem,1), dtype='object')                 ## Element lateral gravitation coefficient
## Element coordinates
coords[0,0] = np.array([[0,0],
                        [20,10],
                        [0,10]], dtype='float') ## Element 1
              
coords[1,0] = np.array([[0,0],
                        [20,0],
                        [20,10]], dtype='float') ## Element 2
## Element connection
connect[0,0] = np.array([0,2,1])
connect[1,0] = np.array([0,3,2])
                   
## Element properties             
t[0,0], t[1,0] = 1, 1
E[0,0], E[1,0] = 30E6, 30E6
nuy[0,0], nuy[1,0] = 0.3, 0.3
gamma[0,0], gamma[1,0] = 100, 100
c[0,0], c[1,0] = 0, 0 
# Load vector initializations
## Body loads
bdloads = np.empty((nelem,1), dtype='object') 
## Distributed loads 
sloads = np.empty((nelem,1), dtype='object') 
## Concentrated loads
cloads = np.empty((nelem,1), dtype='object') 

# Construction of linear natural coordinates
def local_stiffness():
    for i in range(nelem):
        theta[i,0], theta[i,0][:,0] = np.zeros((3,3)), [1,1,1]
        theta[i,0][:,1:] = coords[i,0][:,0:]
        A[i,0] = 0.5 * np.linalg.det(theta[i,0])
        D[i,0] = (E[i,0]/(1-nuy[i,0]**2)) * np.array([[1, nuy[i,0], 0],
                                                     [nuy[i,0], 1, 0],
                                                     [0, 0, (1-nuy[i,0])/2]])
        abg[i,0] = np.matmul(np.linalg.inv(theta[i,0]), eye(3))
        ## Compute general displacement operator
        u[i,0] = np.matmul(Matrix([1,x,y]).T, abg[i,0]) 
        N[i,0],N[i,0][0,0::2],N[i,0][1,1::2] = np.zeros((2,6),dtype='object'),u[i,0],u[i,0]
        ## Computation of shape function
        B[i,0] = np.zeros((3,6))
        for j in range(6):
            B[i,0][0,j],B[i,0][1,j] = diff(N[i,0][0,j],x), diff(N[i,0][1,j],y)
            B[i,0][2,0::2], B[i,0][2,1::2] = B[i,0][1,1::2], B[i,0][0,0::2] 
        ## Construction of element stiffness matrix    
        k[i,0] = t[i,0] * A[i,0] * np.matmul(np.matmul(B[i,0].T, D[i,0]), B[i,0])
    return k, A
    
k = local_stiffness()

def Assembly(k):
    for i in range(nelem):
        K[i,0] = np.zeros((nnode*2,nnode*2))
        
        firstu = connect[i,0][0] * 2
        firstv = firstu + 1
        secondu = connect[i,0][1] * 2
        secondv = secondu + 1
        thirdu = connect[i,0][2] * 2
        thirdv = thirdu + 1
        
        DOF[i,0] = np.array([firstu, firstv, secondu, secondv, thirdu, thirdv],dtype='int')
        
        for m,j in enumerate(DOF[i,0]):
            for n,l in enumerate(DOF[i,0]):
                K[i,0][j,l] = K[i,0][j,l] + k[i,0][m,n]
        
    K_add = np.sum(K)
        
                
    return K_add
K_add = Assembly(k)
'''
def Assembly_force_vector(A,t,gamma,c):
    for i in range(nelem):
        bdloads[i,0] = np.zeros((6,1))
        ## in x-direction
        bdloads[i,1::2] = c[i,1]*1/3*A[i,0]*t[i,0]*gamma[i,0]             
         ## in y-direction
        bdloads[i,0::2] = 1/3*A[i,0]*t[i,0]*gamma[i,0]      
    
    return bdloads

bdloads = Assembly_force_vector(A,t,gamma,c)
'''