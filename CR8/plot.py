from main import coords8N,connect,BCs,coordsNEW
from ini import *
import matplotlib.pyplot as plt
fig, ax = plt.subplots()

for i in range(nelem):
   
    plt.plot(coords8N[i][:2,0],coords8N[i][:2,1],'-b')
    plt.plot(coords8N[i][1:3,0],coords8N[i][1:3,1],'-b')
    plt.plot(coords8N[i][2:4,0],coords8N[i][2:4,1],'-b')
    plt.plot([coords8N[i][0,0],coords8N[i][3,0]],[coords8N[i][0,1],coords8N[i][3,1]],'-b')
    
    plt.plot(coordsNEW[:2,0],coordsNEW[:2,1],'--r')
    plt.plot(coordsNEW[1:3,0],coordsNEW[1:3,1],'--r')
    plt.plot(coordsNEW[2:4,0],coordsNEW[2:4,1],'--r')
    plt.plot([coordsNEW[0,0],coordsNEW[3,0]],[coordsNEW[0,1],coordsNEW[3,1]],'--r')
    
    plt.plot(coords8N[i][:,0],coords8N[i][:,1],'ob')
    plt.plot(coordsNEW[:,0],coordsNEW[:,1],'or')

    NodeLabel = connect[0]
    for m,n in enumerate(NodeLabel):
        ax.annotate(int(n), (coords8N[i][m,0],coords8N[i][m,1]),\
        textcoords='offset points',xytext=(8,8),ha='center',fontsize=12,color='blue')
        fig.patch.set_visible('False')   
'''        
for j in BCs[:,0]:
    if BCs[j,1] == 0:
        plt.plot(coords8N[i][j,0],coords8N[i][j,1],'r>',markersize=10)
    if BCs[j,1] == 1:
        plt.plot(coords8N[i][j,0],coords8N[i][j,1],'r^',markersize=10)    
    print(j)
'''    

#ax.axis('off')        
plt.show()   