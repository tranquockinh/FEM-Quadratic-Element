from sympy import *
import numpy as np
from sympy.interactive.printing import init_printing
init_printing(use_unicode=False, wrap_line=False)

# Coordinate definition and initializations
x,y = symbols('x y')
nelem = 2                                               ## Numnber of elements
elemNode = 4                                            ## Node per element
nnode = 6                                               ## Numnber of nodes
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
Princ_Stress =  np.empty((nelem,1), dtype='object')     ## Element principal stresses 

## Element coordinates
coords[0,0] = np.array([[0,0],
                        [0,10],
                        [10,10],
                        [10,0]], dtype='float') ## Element 1
              
coords[1,0] = np.array([[10,0],
                        [10,10],
                        [20,10],
                        [20,0]], dtype='float') ## Element 2
## Element connection
connect[0,0] = np.array([0,1,2,3])
connect[1,0] = np.array([3,2,5,4])
                   
## Element properties             
t[0,0], t[1,0] = 1, 1
E[0,0], E[1,0] = 30E6, 30E6
nuy[0,0], nuy[1,0] = 0.3, 0.3
gamma[0,0], gamma[1,0] = 0, 0
c[0,0], c[1,0] = 0, 0 

# Load vector initializations
## Distributed loads (equivalent concentrated loads)
sfloads[0,0] = np.array([0,0,0,0,0,0,0,0])[:,np.newaxis]
sfloads[1,0] = np.array([0,0,0,0,5000,0,5000,0])[:,np.newaxis]
print()
## Concentrated loads
cloads[0,0] = np.array([0,0,0,0,0,0,0,0])[:,np.newaxis]
cloads[1,0] = np.array([0,0,0,0,0,0,0,0])[:,np.newaxis]

# Boundary conditions
BCs = np.array([[0, 0, 0],
                [0, 1, 0],
                [1, 0, 0],
                [1, 1, 0]])

# Construction of linear natural coordinates
def local_stiffness():
    for i in range(nelem):
        theta[i,0], theta[i,0][:,0] = np.zeros((elemNode,elemNode)), [1,1,1,1]
        theta[i,0][:,1:3] = coords[i,0][:,0:]
        theta[i,0][:,-1] = coords[i,0][:,0] * coords[i,0][:,1]
        
        
        ## For plane stress sigmaz = 0
        
        D[i,0] = (E[i,0]/(1-nuy[i,0]**2)) * np.array([[1, nuy[i,0], 0],
                                                     [nuy[i,0], 1, 0],
                                                     [0, 0, (1-nuy[i,0])/2]])
        '''
        ## For plane stress deltaz = 0  
        D[i,0] = (E[i,0]/((1+nuy[i,0])*(1-2*nuy[i,0]))) * np.array([[1-nuy[i,0], nuy[i,0], 0],
                                                                   [nuy[i,0], 1-nuy[i,0], 0],
                                                                   [0, 0, (1-2*nuy[i,0])/2]])
        '''                                                          
        abcd[i,0] = np.matmul(np.linalg.inv(theta[i,0]), eye(elemNode))
        ## Compute general displacement operator
        Lu[i,0] = np.matmul(Matrix([1,x,y,x*y]).T, abcd[i,0]) 
        N[i,0],N[i,0][0,0::2],N[i,0][1,1::2] = np.zeros((2,2*elemNode),dtype='object'),Lu[i,0],Lu[i,0]
        ## Computation of shape function
        B[i,0] = np.zeros((3,2*elemNode),dtype='object')
        for j in range(2*elemNode):
            B[i,0][0,j],B[i,0][1,j] = diff(N[i,0][0,j],x), diff(N[i,0][1,j],y)
            B[i,0][2,0::2], B[i,0][2,1::2] = B[i,0][1,1::2], B[i,0][0,0::2] 
        ## Construction of element stiffness matrix    
        k[i,0] = t[i,0] * np.matmul(np.matmul(B[i,0].T, D[i,0]), B[i,0])   
        for m in range(2*elemNode):
            for l in range(2*elemNode):
                xlow = coords[i,0][1,0]
                xup = coords[i,0][2,0] 
                ylow = coords[i,0][0,1]
                yup = coords[i,0][1,1] 
                A[i,0] = (xup-xlow) * (yup-ylow)
                k[i,0][m,l] = integrate(integrate(k[i,0][m,l],(x,xlow,xup)),(y,ylow,yup))
    return k,N,B,D,A
    
k,N,B,D,A = local_stiffness()

def Assembly_Stiffness_Matrix(k):
    for i in range(nelem):
        ## Gobal stiffness matrix 
        K[i,0] = np.zeros((nDOF,nDOF))
        
        firstu = connect[i,0][0] * 2
        firstv = firstu + 1
        secondu = connect[i,0][1] * 2
        secondv = secondu + 1
        thirdu = connect[i,0][2] * 2
        thirdv = thirdu + 1
        fourthu = connect[i,0][3] * 2
        fourthv = fourthu + 1
        DOF[i,0] = np.array([firstu, firstv, secondu, secondv, thirdu, thirdv, fourthu, fourthv])
        
        for m,j in enumerate(DOF[i,0]):
            for n,l in enumerate(DOF[i,0]):
                K[i,0][j,l] = K[i,0][j,l] + k[i,0][m,n]      
    K_add = np.sum(K)
    return K_add
K_add = Assembly_Stiffness_Matrix(k)

def Assembly_force_vector(A,t,gamma,c,sfloads,cloads):
    for i in range(nelem):
        bdloads[i,0] = np.zeros((nDOF,1))
        ## in y-direction
        bdloads[i,0][1::2,0] = 1/3*A[i,0]*t[i,0]*gamma[i,0]      
        ## in x-direction
        bdloads[i,0][0::2,0] = c[i,0]*bdloads[i,0][0::2,0]  
        
        F[i,0] = np.zeros((nnode*2,1))
        
        firstu = connect[i,0][0] * 2
        firstv = firstu + 1
        secondu = connect[i,0][1] * 2
        secondv = secondu + 1
        thirdu = connect[i,0][2] * 2
        thirdv = thirdu + 1
        fourthu = connect[i,0][3] * 2
        fourthv = fourthu + 1
        DOF[i,0] = np.array([firstu, firstv, secondu, secondv, thirdu, thirdv, fourthu, fourthv])
        
        for m,j in enumerate(DOF[i,0]):
            F[i,0][j,0] = F[i,0][j,0] + bdloads[i,0][m,0] + sfloads[i,0][m,0] + cloads[i,0][m,0]      
        ## Compute the assembled load
    F_add = np.sum(F)
    return F_add

F_add = Assembly_force_vector(A,t,gamma,c,sfloads,cloads)
def Solver(K,F):
    U = np.zeros((nDOF, 1))
    nPresDOF = len(BCs[:,0])
    nFreeDOF = nDOF - nPresDOF
    idxPresDOF = np.zeros((nPresDOF, 1), dtype='int')
    idxFreeDOF = np.zeros((nFreeDOF, 1), dtype='int')
     
    for i in range(nPresDOF):
        nodePres = BCs[i,0]
        xORy = BCs[i,1]
        Value = BCs[i,2]
        idxPres = nodePres*2 + xORy
        U[idxPresDOF] = Value
        idxPresDOF[i,0] = idxPres
    
    idxFreeDOF = np.arange(0,nDOF)
    idxFreeDOF = np.delete(idxFreeDOF, idxPresDOF)
    
    K_Free = np.zeros((len(idxFreeDOF),len(idxFreeDOF)))
    u_Free = np.zeros((len(idxFreeDOF),1))
    F_Free = np.zeros((len(idxFreeDOF),1))
    for m,i in enumerate(idxFreeDOF):
        for n,j in enumerate(idxFreeDOF):
            K_Free[m,n] = K[i,j]
            F_Free[m] = F[i]
            
    K_Pres = np.zeros((len(idxPresDOF),len(idxPresDOF)))
    F_Pres = np.zeros((len(idxPresDOF),1))
    u_Pres = np.zeros((len(idxPresDOF),1))
    for m,i in enumerate(idxPresDOF):
        for n,j in enumerate(idxPresDOF):
            K_Pres[m,n] = K[i,j]
            F_Pres[m] = F[i]
            u_Pres[m] = U[i]

    
    K_PresFree = np.zeros((len(idxPresDOF),len(idxFreeDOF)))
    K_FreePres = np.zeros((len(idxFreeDOF),len(idxPresDOF)))
    for m,i in enumerate(idxPresDOF):
        for n,j in enumerate(idxFreeDOF):
            K_PresFree[m,n] = K[i,j]
            K_FreePres[n,m] = K[j,i]

    u_Free = np.matmul(np.linalg.inv(K_Free),(F_Free - np.matmul(K_FreePres,u_Pres)))
    p_Pres = np.matmul(K_PresFree,u_Free) + np.matmul(K_Pres,u_Pres)
    ## COncatenate to full vector of displacement and forces
    for i,j in enumerate(idxFreeDOF):
        U[j] = u_Free[i]
    for i,j in enumerate(idxPresDOF):
        F[j] = p_Pres[i]
    EU = np.zeros((2,nnode))
    EF = np.zeros((2,nnode))
    for i in range(nnode):
        EU[:,i] = U[2*i:2*i+2,0]
        EF[:,i] = F[2*i:2*i+2,0]
    #print('Global displacement is:\n',U)
    print()
    #print('Global force vector is:\n',F)
    print()
    print('Element displacement is:\n',EU.T)
    print()
    print('Element force vector is:\n',EF.T)
    print()
    return U,EU,F,EF,K_Free

U,EU,F,EF,K_Free = Solver(K_add,F_add)   

print(np.linalg.inv(K_Free))