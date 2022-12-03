from ini import *

def Global_Stiffness_Matrix(k,connect,K):
    for i in range(nelem):
        fitstu = connect[0] * 2
        firstv = fitstu + 1
        secondu = connect[2] * 2
        secondv = secondu + 1
        thirdu = connect[4] * 2
        thirdv = thirdu + 1
        fourthu = connect[6] * 2
        fourthv = fourthu + 1
        fifthu = connect[1] * 2
        fifthv = fifthu + 1
        sixthu = connect[3] * 2
        sixthv = sixthu + 1
        seventhu = connect[5] * 2
        seventhv = seventhu + 1
        eighthu = connect[7] * 2
        eighthv = eighthu + 1
        
        DOF = np.array([fitstu,firstv,secondu,secondv,thirdu,thirdv,fourthu,fourthv,\
                            fifthu,fifthv,sixthu,sixthv,seventhu,seventhv,eighthu,eighthv])
        
        
        for l,m in enumerate(DOF):
            for j,n in enumerate(DOF):
               
               K[m,n] = K[m,n] + k[l,j] 
               
        return K
    