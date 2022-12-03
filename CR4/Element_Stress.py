from sympy import *
import numpy as np
from Initial import *

def Element_Stress(D_Store,B_Store,U,connect):

    EStress = np.zeros((3,nelem),dtype='float')
    
    ## Go back local displacement vector
    for i in range(nelem):
    ## Each element has 4 nodes
        u[i] = np.zeros((2*elemNode,1))
        Stress[i] = np.zeros((3,1))
        
        firstu = connect[i][0] * 2
        firstv = firstu + 1
        secondu = connect[i][1] * 2
        secondv = secondu + 1
        thirdu = connect[i][2] * 2
        thirdv = thirdu + 1
        fourthu = connect[i][3] * 2
        fourthv = fourthu + 1
        
        DOF = np.array([firstu, firstv, secondu, secondv, thirdu, thirdv, fourthu, fourthv])    
        
        for j,l in enumerate(DOF):
            u[i][j,0] = U[l] 
        
        Stress[i][:,0] = np.matmul(np.matmul(D_Store[i],B_Store[i]),u[i])
        EStress[:,i] = Stress[i][:,0]
        
    return Stress,EStress
    