from ini import *
from stiffness import CR8_Stiffness

def Global_Force_Vector(connect,fs,fb,Fs,Fb):
    
    firstu = connect[0] * 2
    firstv = firstu + 1
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
    
    DOF = np.array([firstu,firstv,secondu,secondv,thirdu,thirdv,fourthu,fourthv,\
                        fifthu,fifthv,sixthu,sixthv,seventhu,seventhv,eighthu,eighthv])
       
    for j,l in enumerate(DOF):  
        Fs[l] = fs[j]
        Fb[l] = fb[j]
    
    F_global = Fs + fb
    
    return F_global