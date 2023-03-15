import numpy
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.ticker import LinearLocator
from matplotlib import cm
import numpy as np
from matplotlib import pyplot as plt
from numpy.linalg import eig

def Eigenvalue (mass, a, b, front_stiff, rear_stiff, yaw_inertia,):   
    velocity = 1 # in m/s

    #Entries of Matrix A [d e , f g]

    eig_1 = -10
    eig_2 = -10
    

    while True: 
        d = ((-(front_stiff + rear_stiff))/(mass*velocity)) 
        e = ((-(a*front_stiff-b*rear_stiff))/(mass*velocity))-velocity
        f = ((-(a*front_stiff-b*rear_stiff))/(yaw_inertia*velocity))
        g = ((-(a**2*front_stiff+b**2*rear_stiff)/(yaw_inertia*velocity)))
        
        matrix = np.array([[d, e], [f, g]])
        eigenvalues = eig(matrix)
        velocity +=1
        if ((eigenvalues > 0).any()):
            break

    return velocity

velocity = Eigenvalue(1400, 1.14, 1.33, 25000, 21000, 2420)
print(velocity*3.6)



