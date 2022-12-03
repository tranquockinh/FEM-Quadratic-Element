from Initial import *
from EvalE import C4R_Stiffness

coords[0] = np.array([[3,2],
                      [5,2],
                      [5,4],
                      [3,4]], dtype='float')
'''                      
coords[1] = np.array([[3,2],
                      [5,2],
                      [5,4],
                      [3,4]], dtype='float')
'''                      
## Material properties
E[0] = 30E6
nuy[0] = 0.25
thickness[0] = 1

for i in range(nelem):

    k[i] = C4R_Stiffness(coords[i],E[i],nuy[i],thickness[i]) 
       

print(k[i])

