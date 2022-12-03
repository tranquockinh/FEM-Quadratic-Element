from ini import *

def CR8_Stiffness(coords8N,E,nuy,thickness):
           
    N1 = 0.25 * (1-s)*(1-t)*(-s-t-1)
    N2 = 0.25 * (1+s)*(1-t)*(s-t-1)
    N3 = 0.25 * (1+s)*(1+t)*(s+t-1)
    N4 = 0.25 * (1-s)*(1+t)*(-s+t-1)
    # Mid-point
    N5 = 0.5*(1-t)*(1+s)*(1-s)
    N6 = 0.5*(1+t)*(1+s)*(1-t)
    N7 = 0.5*(1+t)*(1+s)*(1-s)
    N8 = 0.5*(1-s)*(1+t)*(1-t)
   
    Ns = np.zeros((2,16),dtype='object')
    Nt = np.zeros((2,16),dtype='object')
    # Shape function matrix
    ShapeFunc = sp.Matrix([[N1,0,N2,0,N3,0,N4,0,N5,0,N6,0,N7,0,N8,0],
                           [0,N1,0,N2,0,N3,0,N4,0,N5,0,N6,0,N7,0,N8]])
    
    for i in range(len(ShapeFunc[:,0])):
        for j in range(len(ShapeFunc[0,:])):
            Ns[i,j] = diff(ShapeFunc[i,j],s)
            Nt[i,j] = diff(ShapeFunc[i,j],t)      
            
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
   
    D1_11 = diff(y[0],t)
    D1_22 = -diff(x[0],t)  
    D2_11 = -diff(y[0],s)
    D2_22 = diff(x[0],s)
      
    Dwrtt = np.array([[D1_11, 0],
                     [0, D1_22],
                     [D1_22, D1_11]])
    Dwrts = np.array([[D2_11, 0],
                     [0, D2_22],
                     [D2_22, D2_11]])
    
    # Compute Jacobian matrix
    JacobianMat = np.zeros((2,2),dtype='object')
    JacobianMat[0,0] = diff(x[0],s)
    JacobianMat[0,1] = diff(y[0],s)
    JacobianMat[1,0] = diff(x[0],t)
    JacobianMat[1,1] = diff(y[0],t)
    JacobianMat = Matrix(JacobianMat)
    detJ = JacobianMat.det()
    
    B = (1/detJ)*(np.matmul(Dwrtt,Ns) + np.matmul(Dwrts,Nt))
    B = Matrix(B)
    
    
    # Compute B matrix at Gausian points
    ## Gaussian points
    sGP = np.array([-1/np.sqrt(3),1/np.sqrt(3)])
    tGP = np.array([-1/np.sqrt(3),1/np.sqrt(3)])
    W = np.array([1,1])
    k = np.zeros((elemNode*2,elemNode*2),dtype='float')
    for m in range(2):
        for n in range(2):
            detJ = detJ.subs([(s,sGP[m]),(t,tGP[n])])
            B = B.subs([(s,sGP[m]),(t,tGP[n])])
            w1 = W[m]
            w2 = W[n]
            k = k + thickness * w1 * w2 * np.matmul(np.matmul(B.T,D),B)
            
    return k,ShapeFunc,detJ,B,D