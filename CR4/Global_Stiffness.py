from Initial import *

def Global_CR4_Stiffness_Matrix(k,connect,K):
    
    fitstu = connect[0] * 2
    firstv = fitstu + 1
    secondu = connect[1] * 2
    secondv = secondu + 1
    thirdu = connect[2] * 2
    thirdv = thirdu + 1
    fourthu = connect[3] * 2
    fourthv = fourthu + 1
    
    DOF = np.array([fitstu,firstv,secondu,secondv,thirdu,thirdv,fourthu,fourthv])
    
    for l,m in enumerate(DOF):
        for j,n in enumerate(DOF):
           
           K[m,n] = K[m,n] + k[l,j] 
           
    return K
    