import numpy
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.ticker import LinearLocator
from matplotlib import cm
import numpy as np
from matplotlib import pyplot as plt



def Euler (iteration, grid, mass, a, b, front_stiff, rear_stiff, yaw_inertia, velocity_x, delta):   
    velocity = velocity_x / 3.6 
    d = ((-(front_stiff + rear_stiff))/(mass*velocity))
    e = ((-(a*front_stiff-b*rear_stiff))/(mass*velocity))-velocity
    f = ((-(a*front_stiff-b*rear_stiff))/(yaw_inertia*velocity))
    g = ((-(a**2*front_stiff+b**2*rear_stiff)/(yaw_inertia*velocity)))

    h = (front_stiff/mass)
    p = ((a*front_stiff)/yaw_inertia)
    
    y = []
    psi = []
    time = []
    lateral_acceleration = []

    y = [0]*iteration
    psi = [0]*iteration
    time = [0]*iteration
    lateral_acceleration = [0]*iteration


    for i in range (1,iteration):
            y[i] = y[i-1] + (grid)*(d*y[i-1]+e*psi[i-1]+h*delta)
            psi[i] = psi[i - 1] + (grid)*(f*y[i-1]+g*psi[i-1]+p*delta)

            time[i] = time[i-1]+ grid
            lateral_acceleration[i] = d*y[i]+e*psi[i]+h*delta

    return y, psi, lateral_acceleration, time
    


def rk (iteration, grid, mass, a, b, front_stiff, rear_stiff, yaw_inertia, car_velocity, delta): 
    velocity = car_velocity / 3.6 
    d = ((-(front_stiff + rear_stiff))/(mass*velocity))
    e = ((-(a*front_stiff-b*rear_stiff))/(mass*velocity))-velocity
    f = ((-(a*front_stiff-b*rear_stiff))/(yaw_inertia*velocity))
    g = ((-(a**2*front_stiff+b**2*rear_stiff)/(yaw_inertia*velocity)))

    h = (front_stiff/mass)
    p = ((a*front_stiff)/yaw_inertia)
    
    y = []
    psi = []
    time = []
    lateral_acceleration = []

    psi_pos = []
    y_pos = []

    velocity_x = []
    velocity_y= [] 

    position_x = []
    position_y = []



    y = [0]*iteration
    psi = [0]*iteration
    time = [0]*iteration
    lateral_acceleration = [0]*iteration
    y_pos = [0]*iteration 
    psi_pos = [0]*iteration

    velocity_x = [0]*iteration
    velocity_y= [0]*iteration 

    position_x = [0]*iteration
    position_y = [0]*iteration




    for i in range (1,iteration):
        f_one_y = d*y[i-1]+e*psi[i-1]+h*delta
        f_two_y = d*(y[i-1]+grid/2)+e*(psi[i-1]+(grid/2)*f_one_y)+h*delta
        f_three_y = d*(y[i-1]+grid/2)+e*(psi[i-1]+(grid/2)*f_two_y)+h*delta
        f_four_y = d*(y[i-1]+grid)+e*(psi[i-1]+grid*f_three_y) + h*delta
        y[i] = y[i-1] + (grid/6)*(f_one_y+2*f_two_y+2*f_three_y+f_four_y)

    
        f_one_psi = f*y[i-1]+g*psi[i-1]+p*delta
        f_two_psi = f*(y[i-1]+grid/2)+g*(psi[i-1]+(grid/2)*f_one_psi)+p*delta
        f_three_psi = f*(y[i-1]+grid/2)+g*(psi[i-1]+(grid/2)*f_two_psi)+p*delta
        f_four_psi = f*(y[i-1]+grid)+g*(psi[i-1]+grid*f_three_psi) + p*delta
        psi[i] = psi[i-1] + (grid/6)*(f_one_psi+2*f_two_psi+2*f_three_psi+f_four_psi)

        time[i] = time[i-1]+ grid
        lateral_acceleration[i] = d*y[i]+e*psi[i]+h*delta

        ## Finding Psi tp plug into position equations
        y_pos[i] = y_pos[i-1] + grid*(y[i-1]) #may need to change indexing since it's calculating y_pos by using the previous y[i-1], even though y[i] is already calculated
        psi_pos[i] = psi_pos[i-1] + grid*(psi[i-1]) #euler's method

        velocity_x [i] = velocity*np.cos(psi_pos[i])-(y[i]+a*psi[i])*np.sin(psi_pos[i]) #velocities calculated using euler's method
        velocity_y [i] = (y[i]+a*psi[i])*np .cos(psi_pos[i]) + velocity*np.sin(psi_pos[i])
        
        position_x [i] = position_x[i-1] + grid*( velocity_x [i-1]) #absolute position calculated using eulers method
        position_y [i] = position_y[i-1] + grid*( velocity_y [i-1])

    


    return y, psi, lateral_acceleration,position_x,position_y, velocity_x, velocity_y, time

y_200, psi_200, lat_accel_200, position_x_200, position_y_200,velocity_x_200, velocity_y_200, time_200 = rk(5000, 0.01, 1400, 1.14, 1.33, 25000, 21000, 2420, 200, 0.1)
y_220, psi_220, lat_accel_220, position_x_220, position_y_220,velocity_x_220, velocity_y_220, time_220 = rk(5000, 0.01, 1400, 1.14, 1.33, 25000, 21000, 2420, 220, 0.1)
y_234, psi_234, lat_accel_234, position_x_234, position_y_234,velocity_x_234, velocity_y_234, time_234 = rk(5000, 0.01, 1400, 1.14, 1.33, 25000, 21000, 2420, 234, 0.1)
y_250, psi_250, lat_accel_250, position_x_250, position_y_250,velocity_x_250, velocity_y_250, time_250 = rk(5000, 0.01, 1400, 1.14, 1.33, 25000, 21000, 2420, 250, 0.1)
y_300, psi_300, lat_accel_300, position_x_300, position_y_300,velocity_x_300, velocity_y_300, time_300 = rk(5000, 0.01, 1400, 1.14, 1.33, 25000, 21000, 2420, 300, 0.1)





fig, axs = plt.subplots(2, figsize=(5,5))
axs[0].set_title('Lateral acceleration vs time')
axs[0].set_xlabel("Time (s)")
axs[0].set_ylabel("Lateral Acceleration (m/s^2)")
axs[0].plot(time_200, lat_accel_200, label = "200 km/h")
axs[0].plot(time_220, lat_accel_220, label = "220 km/h")
axs[0].plot(time_234, lat_accel_234, label = "234 km/h")
axs[0].plot(time_250, lat_accel_250, label = "250 km/h")
axs[0].plot(time_300, lat_accel_300, label = "300 km/h")
axs[0].legend()
axs[0].legend(loc="lower left")

axs[1].set_title('Yaw rate vs time')
axs[1].set_xlabel("Time (s)")
axs[1].set_ylabel("Yaw rate (rad/s)")
axs[1].plot(time_200, psi_200, label = "200 km/h")
axs[1].plot(time_220, psi_220, label = "220 km/h")
axs[1].plot(time_234, psi_234, label = "234 km/h")
axs[1].plot(time_250, psi_250, label = "250 km/h")
axs[1].plot(time_300, psi_300, label = "300 km/h")
axs[1].legend()
axs[1].legend(loc="upper left")

# axs.set_title('Absolute position of Car')
# axs.set_xlabel("X (m)")
# axs.set_ylabel("Y (m)")
# axs.plot(velocity_x_20, velocity_y_20, label = "20 km/h")

# axs.legend()
# axs.legend(loc="upper right")

plt.tight_layout()
plt.show()