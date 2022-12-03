from sympy import *
import numpy as np
import sympy as sp
#from sympy.interactive.printing import init_printing
#init_printing(use_unicode=False, wrap_line=False)
elemNode = 4
nelem = 1 
## Variable
s,t = symbols('s t')

## Initial objects                                       
coords = np.empty((nelem,1),dtype='object')
D = np.empty((nelem,1),dtype='object')
E = np.empty((nelem,1),dtype='object')
nuy = np.empty((nelem,1),dtype='object')
x = np.empty((nelem,1),dtype='object')
y = np.empty((nelem,1),dtype='object')
JacobianMat = np.empty((nelem,1),dtype='object')
detJ = np.empty((nelem,1),dtype='object')
B = np.empty((nelem,1),dtype='object') 
k = np.empty((nelem,1),dtype='object') 
thickness = np.empty((nelem,1),dtype='object') 


coords[0,0] = np.array([[3,2],
                        [5,2],
                        [5,4],
                        [3,4]], dtype='float')
E[0,0] = 30E6
nuy[0,0] = 0.25
thickness[0,0] = 1
                 
## Gausian points                        
tGP = np.array([-1/np.sqrt(3),1/np.sqrt(3)])
sGP = np.array([-1/np.sqrt(3),1/np.sqrt(3)]) 
w1, w2 = 1, 1
## Shape functions
N1 = 0.25 * (1-s) * (1-t)
N2 = 0.25 * (1+s) * (1-t)
N3 = 0.25 * (1+s) * (1+t)
N4 = 0.25 * (1-s) * (1+t)


ShapeFunc = sp.Matrix([[N1, 0, N2, 0, N3, 0, N4, 0],
                       [0, N1, 0, N2, 0, N3, 0, N4]])
# Coefficient matrix of detJ  
detJMat = 1/8 * sp.Matrix([[0, 1-t, t-s, s-1],
                           [t-1, 0, s+1, -s-t],
                           [s-t, -s-1, 0, t+1],
                           [1-s, s+t, -t-1, 0]])              
for i in range(nelem):    
    
    # Compute elastic coefficient matrix
    D[i,0] = (E[i,0]/(1-nuy[i,0]**2)) * np.array([[1, nuy[i,0], 0],
                                            [nuy[i,0], 1, 0],
                                            [0, 0, (1-nuy[i,0])/2]])
    # Compute global coordinates
    x[i,0] = np.matmul(ShapeFunc[0,0::2],coords[0,0][:,0])
    y[i,0] = np.matmul(ShapeFunc[1,1::2],coords[0,0][:,1]) 
    
    # Compute Jacobian matrix
    JacobianMat[i,0] = np.zeros((2,2),dtype='object')
    JacobianMat[i,0][0,0] = diff(x[i,0][0],s)
    JacobianMat[i,0][0,1] = diff(y[i,0][0],s)
    JacobianMat[i,0][1,0] = diff(x[i,0][0],t)
    JacobianMat[i,0][1,1] = diff(y[i,0][0],t)
    
    # Compute B matrix at Gausian points
    B[i,0] = np.zeros((3,8),dtype='object')
    k[i,0] = np.zeros((8,8),dtype='float')
    
    for m in range(2):
        for n in range(2):
            detJMat = detJMat.subs([(s,sGP[m]),(t,tGP[n])])
            detJ[i,0] = np.matmul(np.matmul(coords[i,0][:,0].T,detJMat),coords[i,0][:,1])
            #JacobianMat[i] = sp.Matrix(JacobianMat[i]).subs([(s,sGP[m]),(t,tGP[n])])
            
            a = JacobianMat[i,0][1,1].subs([(s,sGP[m]),(t,tGP[n])])
            b = JacobianMat[i,0][0,1].subs([(s,sGP[m]),(t,tGP[n])])
            c = JacobianMat[i,0][0,0].subs([(s,sGP[m]),(t,tGP[n])])
            d = JacobianMat[i,0][1,0].subs([(s,sGP[m]),(t,tGP[n])])
            
            N1s, N1t = diff(N1,s).subs([(t,tGP[n])]), diff(N1,t).subs([(s,sGP[m])])
            N2s, N2t = diff(N2,s).subs([(t,tGP[n])]), diff(N2,t).subs([(s,sGP[m])])
            N3s, N3t = diff(N3,s).subs([(t,tGP[n])]), diff(N3,t).subs([(s,sGP[m])])
            N4s, N4t = diff(N4,s).subs([(t,tGP[n])]), diff(N4,t).subs([(s,sGP[m])])
            
            B1 = sp.Matrix([[a*N1s-b*N1t,  0],
                            [0,            c*N1t-d*N1s],
                            [c*N1t-d*N1s,  a*N1s-b*N1t]])
                            
            B2 = sp.Matrix([[a*N2s-b*N2t,  0],
                            [0,            c*N2t-d*N2s],
                            [c*N2t-d*N2s,  a*N2s-b*N2t]])
            
            B3 = sp.Matrix([[a*N3s-b*N3t,  0],
                            [0,            c*N3t-d*N3s],
                            [c*N3t-d*N3s,  a*N3s-b*N3t]])  
            
            B4 = sp.Matrix([[a*N4s-b*N4t,  0],
                            [0,            c*N4t-d*N4s],
                            [c*N4t-d*N4s,  a*N4s-b*N4t]])  
                            
            B[i,0][:,:2] = B1
            B[i,0][:,2:4] = B2
            B[i,0][:,4:6] = B3
            B[i,0][:,6:8] = B4
            B[i,0] = (1/detJ) * B[i,0]
            
            k[i,0] = k[i,0] + thickness[i,0]*w1*w2*np.matmul(np.matmul(B[i,0].T,D[i,0]),B[i,0])
            k[i,0] = np.array(k[i,0],dtype='float')

print(np.linalg.inv(k[0,0]))
