from ini import *

def CR8_Stiffness(coords8N,E,nuy,thickness):
    
    N1 = 0.25 * (1-s)*(1-t)*(-s-t-1)
    N2 = 0.25 * (1+s)*(1-t)*(s-t-1)
    N3 = 0.25 * (1+s)*(1+t)*(s+t-1)
    N4 = 0.25 * (1-s)*(1+t)*(-s+t-1)
    N5 = 0.5 * (1-t)*(1+s)*(1-s)
    N6 = 0.5 * (1+t)*(1+s)*(1-t)
    N7 = 0.5 * (1+t)*(1+s)*(1-s)
    N8 = 0.5 * (1-s)*(1+t)*(1-t)
   
    Ns = np.zeros((2,16),dtype='object')
    Nt = np.zeros((2,16),dtype='object')
    
    # Shape function matrix
    ShapeFunc = sp.Matrix([[N1,0,N2,0,N3,0,N4,0,N5,0,N6,0,N7,0,N8,0],
                           [0,N1,0,N2,0,N3,0,N4,0,N5,0,N6,0,N7,0,N8]])
    
     
            
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
    x = np.matmul(ShapeFunc[0,0::2],coords8N[:,0])
    y = np.matmul(ShapeFunc[1,1::2],coords8N[:,1])   
    # Compute Jacobian matrix
    JacobianMat = np.zeros((2,2),dtype='object')
    JacobianMat[0,0] = diff(x[0],s)
    JacobianMat[0,1] = diff(y[0],s)
    JacobianMat[1,0] = diff(x[0],t)
    JacobianMat[1,1] = diff(y[0],t)
    JacobianMat = Matrix(JacobianMat)
    detJ = JacobianMat.det()
    '''
    D1_11 = diff(y[0],t)
    D1_22 = -diff(x[0],t)  
    D2_11 = -diff(y[0],s)
    D2_22 = diff(x[0],s)
    
    for i in range(len(ShapeFunc[:,0])):
        for j in range(len(ShapeFunc[0,:])):
            Ns[i,j] = diff(ShapeFunc[i,j],s)
            Nt[i,j] = diff(ShapeFunc[i,j],t)     
            
    Dwrtt = np.array([[D1_11, 0],
                     [0, D1_22],
                     [D1_22, D1_11]])
                     
    Dwrts = np.array([[D2_11, 0],
                     [0, D2_22],
                     [D2_22, D2_11]])
    
    B = np.matmul(Dwrtt,Ns) + np.matmul(Dwrts,Nt)
    B = Matrix(B)
    '''
    
    B = np.zeros((3,16),dtype='object')
    k = np.zeros((8,16),dtype='float')
    # Compute B matrix at Gausian points
    ## Gaussian points
    sGP = np.array([-1/np.sqrt(3),1/np.sqrt(3)])
    tGP = np.array([-1/np.sqrt(3),1/np.sqrt(3)])
    W = np.array([1,1])
    k = np.zeros((elemNode*2,elemNode*2),dtype='float')
    for m in range(2):
        for n in range(2):
            detJ = detJ.subs([(s,sGP[m]),(t,tGP[n])])
            '''
            B =  B.subs([(s,sGP[m]),(t,tGP[n])])
            B = (1/detJ) * B
            '''
            a = JacobianMat[1,1].subs([(s,sGP[m]),(t,tGP[n])])
            b = JacobianMat[0,1].subs([(s,sGP[m]),(t,tGP[n])])
            c = JacobianMat[0,0].subs([(s,sGP[m]),(t,tGP[n])])
            d = JacobianMat[1,0].subs([(s,sGP[m]),(t,tGP[n])])
            
            N1s, N1t = diff(N1,s).subs([(t,tGP[n])]), diff(N1,t).subs([(s,sGP[m])])
            N2s, N2t = diff(N2,s).subs([(t,tGP[n])]), diff(N2,t).subs([(s,sGP[m])])
            N3s, N3t = diff(N3,s).subs([(t,tGP[n])]), diff(N3,t).subs([(s,sGP[m])])
            N4s, N4t = diff(N4,s).subs([(t,tGP[n])]), diff(N4,t).subs([(s,sGP[m])])
            N5s, N5t = diff(N5,s).subs([(t,tGP[n])]), diff(N5,t).subs([(s,sGP[m])])
            N6s, N6t = diff(N6,s).subs([(t,tGP[n])]), diff(N6,t).subs([(s,sGP[m])])
            N7s, N7t = diff(N7,s).subs([(t,tGP[n])]), diff(N7,t).subs([(s,sGP[m])])
            N8s, N8t = diff(N8,s).subs([(t,tGP[n])]), diff(N8,t).subs([(s,sGP[m])])
            
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
            
            B5 = sp.Matrix([[a*N5s-b*N5t,  0],
                            [0,            c*N5t-d*N5s],
                            [c*N5t-d*N5s,  a*N5s-b*N5t]])
                            
            B6 = sp.Matrix([[a*N6s-b*N6t,  0],
                            [0,            c*N6t-d*N6s],
                            [c*N6t-d*N6s,  a*N6s-b*N6t]])
            
            B7 = sp.Matrix([[a*N7s-b*N7t,  0],
                            [0,            c*N7t-d*N7s],
                            [c*N7t-d*N7s,  a*N7s-b*N7t]])  
            
            B8 = sp.Matrix([[a*N8s-b*N8t,  0],
                            [0,            c*N8t-d*N8s],
                            [c*N8t-d*N8s,  a*N8s-b*N8t]])              
            
            B[:,:2] = B1
            B[:,2:4] = B2
            B[:,4:6] = B3
            B[:,6:8] = B4
            B[:,8:10] = B5
            B[:,10:12] = B6
            B[:,12:14] = B7
            B[:,14:16] = B8
            
            B = Matrix(B).subs([(s,sGP[m]),(t,tGP[n])])
            B = (1/detJ * B)
            w1 = W[m]
            w2 = W[n]
            
            k = k + thickness * w1 * w2 * np.matmul(np.matmul(B.T,D),B) * detJ 
            
    return k,ShapeFunc,detJ,B,D,N1s