import numpy
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.ticker import LinearLocator
from matplotlib import cm
import numpy as np
from matplotlib import pyplot as plt
from numpy.linalg import eig

def Eigenvalue (mass, a, b, front_stiff, rear_stiff, yaw_inertia,):   
    velocity = 3000 # in m/s

    #Entries of Matrix A [d e , f g]

    eig_1 = -10
    eig_2 = -10
    

    while True: 
        #Matrix entries
        d = ((-(front_stiff + rear_stiff))/(mass*velocity)) 
        e = ((-(a*front_stiff-b*rear_stiff))/(mass*velocity))-velocity
        f = ((-(a*front_stiff-b*rear_stiff))/(yaw_inertia*velocity))
        g = ((-(a**2*front_stiff+b**2*rear_stiff)/(yaw_inertia*velocity)))
        
        matrix = np.array([[d, e], [f, g]]) #Establishing matrix for numpy
        eigenvalues, eigenvectors = eig(matrix) #Calculating the Eigenvalues
        

        
        if ((eigenvalues.real > -0.01).any()): #if one or both of the eigenvalues are positive, stop the loop and return the velocity corresponding with the previous velocity (to ensure that max speed matches 2 negative eigenvalues)
            break
        velocity += 0.1
        print(eigenvalues.real)
       # print(eigenvalues.imag)

    return velocity #Return the max velocity

velocity = Eigenvalue(1400, 1.14, 1.33, 20000, 20000, 2420)
print(velocity*3.6) #convert max velocity in m/s to km/h



