from ini import *
from stiffness import CR8_Stiffness
from Stiffness_Assemby import Global_Stiffness_Matrix
from Force_Assembly import Global_Force_Vector
from Solver import Solver

'''
====================================================
IF: LINEAR POINT: CHANGE COORDINATE POINTS AT CORNER
IF: CURVATURE: CHNAGE ALL COORDINATE POINTS
====================================================
'''

coords[0] = np.array([[3,2],
                      [5,2],
                      [5,4],
                      [3,4]], dtype='float')
                      
# Element connectivity
connect[0] = np.array([0,1,2,3,4,5,6,7])
# Material properties
E[0] = 30E6
nuy[0] = 0.25
thickness[0] = 1
gamma[0] = 100
Ka[0] = 0

# Load
## Suppose everything is zero
for i in range(nelem):
    Ts[i] = np.zeros((elemNode*2,1),dtype='float')
## Then assisgn where are loaded


fs[0] = np.array([0,0,1000,0,500,0,0,0,0,0,0,0,0,0,0,0])[:,np.newaxis]
fb[0] = np.array([0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0])[:,np.newaxis]
# BCs

BCs = np.array([[0,0,0],
                [0,1,0],
                [1,1,0],
                [3,0,0],
                [4,1,0],
                [7,0,0]])
                
    
for i in range(nelem):
    
    coords8N[i] = np.zeros((8,2))
    coords8N[i][:4,:] = coords[i]
    coords8N[i][4,:] = 0.5*(coords[i][0,:] + coords[i][1,:])
    coords8N[i][5,:] = 0.5*(coords[i][1,:] + coords[i][2,:])
    coords8N[i][6,:] = 0.5*(coords[i][2,:] + coords[i][3,:])
    coords8N[i][7,:] = 0.5*(coords[i][3,:] + coords[i][0,:])
    
    k[i],ShapeFunc,detJ,B,D,N1s = CR8_Stiffness(coords8N[i],E[i],nuy[i],thickness[i])
    B_Store[i] = B
    D_Store[i] = D
    ## Compute global stiffness matrix
    K = Global_Stiffness_Matrix(k[i],connect[i],K)
    
    
    '''
    ## Compute body force vector
    UnitWeight[i] = np.zeros((2,1),dtype='object')
    fb[i] = np.zeros((elemNode*2,1),dtype='object')
    UnitWeight[i] = np.array([gamma[i],Ka[i]*gamma[i]])[:,np.newaxis]
    BodyForce = np.matmul(ShapeFunc.T,UnitWeight[i]) * thickness[i] * detJ
    for j in range(elemNode):
        fb[i][j] = integrate(integrate(BodyForce[j,0],(s,-1,1)),(t,-1,1))
    print(fb)    
    ## Compute Traction
    for j in range(StressNode):
        Traction = Ts[i][j*2:j*2+2,0][:,np.newaxis]
        h = thickness[i]
        Edge0 = np.sqrt(np.sum((coords[i][1,:]-coords[i][0,:])**2))
        Edge1 = np.sqrt(np.sum((coords[i][2,:]-coords[i][1,:])**2))
        Edge2 = np.sqrt(np.sum((coords[i][3,:]-coords[i][2,:])**2))
        Edge3 = np.sqrt(np.sum((coords[i][0,:]-coords[i][3,:])**2))
        Length = np.array([Edge0,Edge1,Edge2,Edge3])
        ## initialize local force on each edge
        Edge = np.zeros((6,1))

        if j==0:
            connectN = np.array([0,4,1])
            N = Matrix([[N1,0,N5,0,N2,0],
                        [0,N1,0,N5,0,N2]]).subs([(t,-1)])           
            insideInt = np.matmul(N.T,Traction) * h * 0.5*Length[j]                    
            for l in range(6):               
                Edge[l] = integrate(insideInt[l][0],(s,-1,1))               
        elif j==1:
            connectN = np.array([1,5,2])
            N = Matrix([[N2,0,N6,0,N3,0],
                        [0,N2,0,N6,0,N3,0]]).subs([(s,1)])                        
            insideInt = (np.matmul(N.T,Traction) * h * 0.5*Length[j])
            for l in range(6):                           
                Edge[l] = integrate(insideInt[l][0],(t,-1,1))                 
        elif j==2:           
            connectN = np.array([2,6,3])
            N = Matrix([[N3,0,N7,0,N4,0],
                        [0,N3,0,N7,0,N4]]).subs([(t,1)])                       
            insideInt = (np.matmul(N.T,Traction) * h * 0.5*Length[j])           
            for l in range(6):               
               Edge[l] = integrate(insideInt[l][0],(s,-1,1))
        else:            
            connectN = np.array([3,7,0])
            N = Matrix([[N4,0,N8,0,N1,0],
                        [0,N4,0,N8,0,N1]]).subs([(s,-1)])                       
            insideInt = (np.matmul(N.T,Traction) * h * 0.5*Length[j])
            for l in range(6):               
                Edge[l] = integrate(insideInt[l][0],(t,-1,1)) 
                
        
        Edge = Matrix(Edge)
        '''
        
        
        
    ## Global stiffness matrix assemble
    F_global = Global_Force_Vector(connect[i],fs[i],fb[i],Fs,Fb)    
    
print(F_global)   
U,F = Solver(F_global,K,BCs)
print('NODAL DISPLACEMENT\n',U)
print()
print('NODAL FORCE\n',F)
print()

# Displacement visualization
coordsNEW = coords8N[0] + U.reshape(8,2) * 300
print(coordsNEW)

