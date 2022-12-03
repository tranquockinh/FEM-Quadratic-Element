from Initial import *
from EvalE import C4R_Stiffness
from Global_Stiffness import Global_CR4_Stiffness_Matrix
from EvalF import Global_Force_Vector
from Solver import Solver
from Element_Stress import Element_Stress

# Coordinates
coords[0] = np.array([[3,2],
                      [5,2],
                      [5,4],
                      [3,4]], dtype='float')
coords[1] = np.array([[5,2],
                      [7,2],
                      [7,4],
                      [5,4]], dtype='float')      

# Shape functions
N1 = 0.25 * (1-s) * (1-t)
N2 = 0.25 * (1+s) * (1-t)
N3 = 0.25 * (1+s) * (1+t)
N4 = 0.25 * (1-s) * (1+t) 
                     
# Element connectivity                      
connect[0] = np.array([0,1,2,3],dtype='int')           
connect[1] = np.array([1,4,5,2],dtype='int')        
    
# Material properties
E[0],E[1] = 30E6,30E6 
nuy[0],nuy[1] = 0.25,0.25
thickness[0],thickness[1] = 1, 1
gamma[0], gamma[1] = 0, 0
Ka[0], Ka[1] = 0, 0

# Loads
for i in range(nelem):
    Ts[i] = np.zeros((8,1),dtype='float')
    
Ts[0] = np.array([0,0,0,0,0,0,0,0])[:,np.newaxis]
Ts[1] = np.array([0,0,-2000,0,0,0,0,0])[:,np.newaxis]

# BCs
BCs = np.array([[0,0,0],
                [0,1,0],
                [1,1,0],
                [3,0,0],
                [3,1,0]])
                
for i in range(nelem):
    ## Compute element stiffness matrix
    k[i],ShapeFunc,detJ,B,D = C4R_Stiffness(coords[i],E[i],nuy[i],thickness[i])
    B_Store[i] = B
    D_Store = D
    ## Compute global stiffness matrix 
    K = Global_CR4_Stiffness_Matrix(k[i],connect[i],K)
    ## Compute global force vector
    fs[i] = np.zeros((8,1),dtype='object')
    ## Compute body force vector
    UnitWeight[i] = np.zeros((2,1),dtype='float')
    fb[i] = np.zeros((8,1),dtype='object')
    UnitWeight[i] = np.array([gamma[i],Ka[i]*gamma[i]])[:,np.newaxis]
    BodyForce = np.matmul(ShapeFunc.T,UnitWeight[i]) * thickness[i] * detJ
    for j in range(8):
        fb[i][j] = integrate(integrate(BodyForce[j,0],(s,-1,1)),(t,-1,1))
            
    ## Compute traction        
    for j in range(elemEdge):
        Traction = Ts[i][j*2:j*2+2,0][:,np.newaxis]
        h = thickness[i]
        Edge0 = np.sqrt(np.sum((coords[i][1,:]-coords[i][0,:])**2))
        Edge1 = np.sqrt(np.sum((coords[i][2,:]-coords[i][1,:])**2))
        Edge2 = np.sqrt(np.sum((coords[i][3,:]-coords[i][2,:])**2))
        Edge3 = np.sqrt(np.sum((coords[i][0,:]-coords[i][3,:])**2))
        Length = np.array([Edge0,Edge1,Edge2,Edge3])
        
        ## initialize local force on each edge
        fs0 = np.zeros((8,1))
        fs1 = np.zeros((8,1))
        fs2 = np.zeros((8,1))
        fs3 = np.zeros((8,1))      
        if j==0:
            N = Matrix([[N1,0,N2,0],
                        [0,N1,0,N2]]).subs([(t,-1)])           
            insideInt = np.matmul(N.T,Traction) * h * 0.5*Length[j]                    
            for l in range(4):               
                fs0[:4,0][l] = integrate(insideInt[l][0],(s,-1,1))               
        elif j==1:            
            N = Matrix([[N2,0,N3,0],
                        [0,N2,0,N3]]).subs([(s,1)])                        
            insideInt = (np.matmul(N.T,Traction) * h * 0.5*Length[j])
            for l in range(4):                           
                fs1[2:6,0][l] = integrate(insideInt[l][0],(t,-1,1))                 
        elif j==2:           
            N = Matrix([[N3,0,N4,0],
                        [0,N3,0,N4]]).subs([(t,1)])                       
            insideInt = (np.matmul(N.T,Traction) * h * 0.5*Length[j])           
            for l in range(4):               
                fs2[4:8,0][l] = integrate(insideInt[l][0],(s,-1,1))
        else:            
            N = Matrix([[N4,0,N1,0],
                        [0,N4,0,N1]]).subs([(s,-1)])                       
            insideInt = (np.matmul(N.T,Traction) * h * 0.5*Length[j])
            for l in range(2):               
                fs3[6:8][l] = integrate(insideInt[:2][l][0],(t,-1,1)) 
                fs3[:2][l] = integrate(insideInt[2:][l][0],(t,-1,1))        
        fs[i] = fs0 + fs1 + fs2 + fs3            
        
        ## Assemble force vector
        F_global = Global_Force_Vector(connect[i],fs[i],fb[i],Fs,Fb)    
    
    
U,F = Solver(F_global,K,BCs)
print('NODAL DISPLACEMENT\n',U)
print()
print('NODAL FORCE\n',F)







