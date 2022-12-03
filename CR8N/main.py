from ini import *
from stiffness import CR8_Stiffness
from Stiffness_Assemby import Global_Stiffness_Matrix
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
BCs = np.array([[0,0,0],
                [0,1,0],
                [4,1,0],
                [1,1,0],
                [3,0,0],
                [3,1,0],
                [7,0,0],
                [7,1,0]])   
                
# Element connectivity
connect[0] = np.array([0,4,1,5,2,6,3,7])
# Material properties
E[0] = 30E6
nuy[0] = 0.25
thickness[0] = 1
gamma[0] = 0
Ka[0] = 0

for i in range(nelem):
    
    coords8N[i] = np.zeros((8,2))
    
    coords8N[i][:4,:] = coords[i]
    coords8N[i][4,:] = 0.5*(coords[i][0,:] + coords[i][1,:])
    coords8N[i][5,:] = 0.5*(coords[i][1,:] + coords[i][2,:])
    coords8N[i][6,:] = 0.5*(coords[i][2,:] + coords[i][3,:])
    coords8N[i][7,:] = 0.5*(coords[i][3,:] + coords[i][0,:])
    
    k[i],ShapeFunc,detJ,B,D = CR8_Stiffness(coords8N[i],E[i],nuy[i],thickness[i])
    B_Store[i] = B
    D_Store[i] = D
    
    ## Compute global stiffness matrix
    K = Global_Stiffness_Matrix(k[i],connect[i],K)

K_F = Solver(K, BCs)
print(K_F)

print(np.linalg.inv(K_F))
