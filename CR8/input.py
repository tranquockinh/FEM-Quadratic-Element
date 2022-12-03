from ini import *
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
gamma[0] = 0
Ka[0] = 0

# Load
## Suppose everything is zero
for i in range(nelem):
    Ts[i] = np.zeros((elemNode,1),dtype='float')
## Then assisgn where are loaded
Ts[0] = np.array([0,0,200,0,0,0,0,0])[:,np.newaxis]

# BCs

BCs = np.array([[0,0,0],
                [0,1,0],
                [1,1,0],
                [3,0,0],
                [3,1,0],
                [4,1,0],
                [7,0,0],
                [7,1,0]])