from Initial import *

def Solver(F, K, BCs):
    
    nPresDOF = len(BCs[:,0])
    nFreeDOF = nDOF - nPresDOF
    idxPresDOF = np.zeros((nPresDOF),dtype='int')
    idxFreeDOF =  np.zeros((nFreeDOF),dtype='int')
    
    ## Locate which nodes have prescribed displacement
    for i in range(nPresDOF):
        NodePres = BCs[i,0]
        xORy = BCs[i,1]
        Value = BCs[i,2]
        AtNode = NodePres * 2 + xORy
        
        idxPresDOF[i] = AtNode
        U[AtNode] = Value
        
    idxFreeDOF  = np.arange(0,nDOF)
    idxFreeDOF = np.delete(idxFreeDOF,idxPresDOF)
    
    K_F = np.zeros((nFreeDOF,nFreeDOF))
    F_F = np.zeros((nFreeDOF,1))
    u_F = np.zeros((nFreeDOF,1))
    for m,i in enumerate(idxFreeDOF):
        for n,j in enumerate(idxFreeDOF):
            K_F[m,n] = K[i,j]
            F_F[m] = F[i]
    K_P = np.zeros((nPresDOF,nPresDOF))
    u_P = np.zeros((nPresDOF,1))
    F_P = np.zeros((nPresDOF,1))
    for m,i in enumerate(idxPresDOF):
        for n,j in enumerate(idxPresDOF):
            K_P[m,n] = K[i,j]
            F_P[m] = F[i]
            u_P[m] = U[i]
            
    K_FP = np.zeros((nFreeDOF,nPresDOF))
    K_PF = np.zeros((nPresDOF,nFreeDOF))    
    for m,i in enumerate(idxFreeDOF):
        for n,j in enumerate(idxPresDOF):
            K_FP[m,n] = K[i,j]
            K_PF[n,m] = K[j,i]
    u_F = np.matmul(np.linalg.inv(K_F),(F_F-np.matmul(K_FP,u_P)))
    p_P = np.matmul(K_PF,u_F) + np.matmul(K_P,u_P)
    
    for i,j in enumerate(idxFreeDOF):
        U[j] = u_F[i]
    for i,j in enumerate(idxPresDOF):
        F[j] = p_P[i]
    
    
    return U,F