from Initial import *
from EvalE import C4R_Stiffness

def Global_Force_Vector(connect,fs,fb,Fs,Fb):
    
    fitstu = connect[0] * 2
    firstv = fitstu + 1
    secondu = connect[1] * 2
    secondv = secondu + 1
    thirdu = connect[2] * 2
    thirdv = thirdu + 1
    fourthu = connect[3] * 2
    fourthv = fourthu + 1
    DOF = np.array([fitstu,firstv,secondu,secondv,thirdu,thirdv,fourthu,fourthv])
    
    for j,l in enumerate(DOF):
        Fs[l] = Fs[l] + fs[j]
        Fb[l] = Fb[l] + fb[j]
    F_global = Fs + Fb 
    return F_global