from ini import *
from stiffness import *
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
    Ts[i] = np.zeros((elemNode,1),dtype='float')
## Then assisgn where are loaded
#Ts[0] = np.zeros((16,2))
Ts[0] = np.array([0,0,0,0,2000,0,2000,0,0,0,0,0,0,0,0,0])[:,np.newaxis]

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
        '''
        Edge40 = np.sqrt(np.sum((coords8N[i][4,:]-coords8N[i][0,:])**2))
        Edge14 = np.sqrt(np.sum((coords8N[i][1,:]-coords8N[i][4,:])**2))
        Edge51 = np.sqrt(np.sum((coords8N[i][5,:]-coords8N[i][1,:])**2))
        Edge25 = np.sqrt(np.sum((coords8N[i][2,:]-coords8N[i][5,:])**2))
        Edge62 = np.sqrt(np.sum((coords8N[i][6,:]-coords8N[i][2,:])**2))
        Edge36 = np.sqrt(np.sum((coords8N[i][3,:]-coords8N[i][6,:])**2))
        Edge73 = np.sqrt(np.sum((coords8N[i][7,:]-coords8N[i][3,:])**2))
        Edge07 = np.sqrt(np.sum((coords8N[i][0,:]-coords8N[i][7,:])**2))
        Length = np.array([Edge40,Edge14,Edge51,Edge25,Edge62,Edge36,Edge73,Edge07])
        '''
        Edge0 = np.sqrt(np.sum((coords[i][1,:]-coords[i][0,:])**2))
        Edge1 = np.sqrt(np.sum((coords[i][2,:]-coords[i][1,:])**2))
        Edge2 = np.sqrt(np.sum((coords[i][3,:]-coords[i][2,:])**2))
        Edge3 = np.sqrt(np.sum((coords[i][0,:]-coords[i][3,:])**2))
        Length = np.array([Edge0,Edge1,Edge2,Edge3])
        ## initialize local force on each edge
        Edge = np.zeros((6,1))
        #Node1 = np.zeros((6,1))
        #Node2 = np.zeros((6,1))
        #Node3 = np.zeros((6,1))   
       
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
                
        
        '''
        Node0 = np.zeros((16,1))
        Node1 = np.zeros((16,1))
        Node2 = np.zeros((16,1))
        Node3 = np.zeros((16,1))   
        Node4 = np.zeros((16,1))
        Node5 = np.zeros((16,1))
        Node6 = np.zeros((16,1))
        Node7 = np.zeros((16,1))   
        
        N04 = np.append(ShapeFunc[:,0:2],ShapeFunc[:,8:10],axis=1)      
        N41 = np.append(ShapeFunc[:,8:10],ShapeFunc[:,2:4],axis=1)      
        N15 = np.append(ShapeFunc[:,2:4],ShapeFunc[:,10:12],axis=1)       
        N52 = np.append(ShapeFunc[:,10:12],ShapeFunc[:,4:6],axis=1)      
        N26 = np.append(ShapeFunc[:,4:6],ShapeFunc[:,12:14],axis=1)      
        N63 = np.append(ShapeFunc[:,12:14],ShapeFunc[:,6:8],axis=1)     
        N37 = np.append(ShapeFunc[:,6:8],ShapeFunc[:,14:16],axis=1)     
        N70 = np.append(ShapeFunc[:,14:16],ShapeFunc[:,0:2],axis=1)
        
        N04 = Matrix(N04).subs([(t,-1)])
        N41 = Matrix(N41).subs([(t,-1)])
        N15 = Matrix(N15).subs([(s,1)])
        N52 = Matrix(N52).subs([(s,1)])
        N26 = Matrix(N26).subs([(t,1)])
        N63 = Matrix(N63).subs([(t,1)])
        N37 = Matrix(N37).subs([(s,-1)])
        N70 = Matrix(N70).subs([(s,-1)])
        
        if j==0:      
            connect[i] = np.array([0,4,1,5,2,6,3,7])
            insideInt = np.matmul(N04.T,Traction) * h * 0.5*Length[j]
            for l in range(4):               
                Node0[:4,0][l] = integrate(insideInt[l][0],(s,-1,0))               
        elif j==1:
            connect[i] = np.array([0,4,1,5,2,6,3,7])
            insideInt = (np.matmul(N41.T,Traction) * h * 0.5*Length[j])     
            for l in range(4):                           
                Node1[2:6,0][l] = integrate(insideInt[l][0],(s,0,1))                 
        elif j==2:
            connect[i] = np.array([0,4,1,5,2,6,3,7])
            insideInt = (np.matmul(N15.T,Traction) * h * 0.5*Length[j])             
            for l in range(4):               
                Node2[4:8,0][l] = integrate(insideInt[l][0],(t,1,0))
        elif j==3:
            connect[i] = np.array([0,4,1,5,2,6,3,7])
            insideInt = (np.matmul(N52.T,Traction) * h * 0.5*Length[j])        
            for l in range(4):               
                Node3[4:8,0][l] = integrate(insideInt[l][0],(t,-1,0))
        elif j==4:
            connect[i] = np.array([0,4,1,5,2,6,3,7])
            insideInt = (np.matmul(N26.T,Traction) * h * 0.5*Length[j])            
            for l in range(4):               
                Node4[6:10,0][l] = integrate(insideInt[l][0],(s,0,1))
        elif j==5:
            connect[i] = np.array([0,4,1,5,2,6,3,7])
            insideInt = (np.matmul(N63.T,Traction) * h * 0.5*Length[j])            
            for l in range(4):               
                Node5[8:12,0][l] = integrate(insideInt[l][0],(s,-1,0))
        elif j==6:
            connect[i] = np.array([0,4,1,5,2,6,3,7])
            insideInt = (np.matmul(N37.T,Traction) * h * 0.5*Length[j])      
            for l in range(4):               
                Node6[12:16,0][l] = integrate(insideInt[l][0],(t,0,1)) 
        else:
            connect[i] = np.array([0,4,1,5,2,6,3,7])
            insideInt = (np.matmul(N70.T,Traction) * h * 0.5*Length[j])            
            print(insideInt)
            for l in range(2):               
                Node7[14:16,0][l] = integrate(insideInt[:2][l][0],(t,-1,0))
                Node7[0:2,0][l] = integrate(insideInt[2:][l][0],(t,-1,0))
        
        '''           
        ## Global stiffness matrix assemble
        F_global = Global_Force_Vector(connectN,fs[i],fb[i],Fs,Fb)    
    
print(F_global)   
U,F = Solver(F_global,K,BCs)
print('NODAL DISPLACEMENT\n',U)
print()
print('NODAL FORCE\n',F)
print()

# Displacement visualization
coordsNEW = coords8N[0] + U.reshape(8,2) * 100
print(coordsNEW)

