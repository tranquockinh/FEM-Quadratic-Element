from ini import *

def Solver(K, BCs):
    
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
   
    
    for m,i in enumerate(idxFreeDOF):
        for n,j in enumerate(idxFreeDOF):
            K_F[m,n] = K[i,j]
           
   
          
   
    return K_F