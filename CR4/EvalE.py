from Initial import *

def C4R_Stiffness(coords,E,nuy,thickness):                          
    ## Gausian points                        
    sGP = np.array([-1/np.sqrt(3),1/np.sqrt(3)]) 
    tGP = np.array([-1/np.sqrt(3),1/np.sqrt(3)])
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
    
        
    # Compute elastic coefficient matrix
    D = (E/(1-nuy**2)) * np.array([[1, nuy, 0],
                                   [nuy, 1, 0],
                                   [0, 0, (1-nuy)/2]])
    # Compute global coordinates
    x = np.matmul(ShapeFunc[0,0::2],coords[:,0])
    y = np.matmul(ShapeFunc[1,1::2],coords[:,1]) 
    
    # Compute Jacobian matrix
    JacobianMat = np.zeros((2,2),dtype='object')
    JacobianMat[0,0] = diff(x[0],s)
    JacobianMat[0,1] = diff(y[0],s)
    JacobianMat[1,0] = diff(x[0],t)
    JacobianMat[1,1] = diff(y[0],t)
    JacobianMat = Matrix(JacobianMat)
    detJ = JacobianMat.det()
    # Compute B matrix at Gausian points
    B = np.zeros((3,8),dtype='object')
    k = np.zeros((8,8),dtype='float')
    
    
    for m in range(2):
        for n in range(2):
            #detJMat = detJMat.subs([(s,sGP[m]),(t,tGP[n])])
            #detJ = np.matmul(np.matmul(coords[:,0].T,detJMat),coords[:,1])
            detJ = detJ.subs([(s,sGP[m]),(t,tGP[n])])
            a = JacobianMat[1,1].subs([(s,sGP[m]),(t,tGP[n])])
            b = JacobianMat[0,1].subs([(s,sGP[m]),(t,tGP[n])])
            c = JacobianMat[0,0].subs([(s,sGP[m]),(t,tGP[n])])
            d = JacobianMat[1,0].subs([(s,sGP[m]),(t,tGP[n])])
            
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
                            
            B[:,:2] = B1
            B[:,2:4] = B2
            B[:,4:6] = B3
            B[:,6:8] = B4
            B = (1/detJ * B)
            
            k = k + thickness*w1*w2*np.matmul(np.matmul(B.T,D),B)*detJ
            
    return k,ShapeFunc,detJ,B,D