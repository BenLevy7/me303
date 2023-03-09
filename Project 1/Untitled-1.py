import numpy
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.ticker import LinearLocator
from matplotlib import cm
import numpy as np
from matplotlib import pyplot as plt
grid = 0.01
iteration  = 10000
initial_y = 0
initial_psi = 0

m = 1400
a = 1.14
b = 1.33
front_stiffness = 25000
rear_stiffness = 21000
yaw_inertia = 2420
velocity_x = 200/3.6
delta = 0.1

d = ((-(front_stiffness + rear_stiffness))/(m*velocity_x))
e = ((-(a*front_stiffness-b*rear_stiffness))/(m*velocity_x))-velocity_x
f = ((-(a*front_stiffness-b*rear_stiffness))/(yaw_inertia*velocity_x))
g = ((-(a**2*front_stiffness+b**2*rear_stiffness)/(yaw_inertia*velocity_x)))


h = (front_stiffness/m)
p = ((a*front_stiffness)/yaw_inertia)

time = []
y = []
lateral_acceleration =[]
psi = []
real_solution_y = []
real_solution_psi = []

y_rk=[]
psi_rk = []



time = [0] * iteration
y = [0] * iteration
psi = [0] * iteration
print(g)

y_rk=[0] * iteration
psi_rk = [0] * iteration
real_solution_y= [0] * iteration
real_solution_psi= [0] * iteration
lateral_acceleration =[0]*iteration



for i in range(1, iteration):


    #Euler's Method
    y[i] = y[i - 1] + (grid)*(d*y[i-1]+e*psi[i-1]+h*delta)
    psi[i] = psi[i - 1] + (grid)*(f*y[i-1]+g*psi[i-1]+p*delta)
    time[i] = time[i - 1] + grid

    #Real Solution
    real_solution_y[i] = -13.0964*numpy.exp(-1.9745*time[i])+24.4684*numpy.exp(-0.9839*time[i])-11.3720
    real_solution_psi[i] = -0.2946*numpy.exp(-1.9745*time[i])-0.6962*numpy.exp(-0.9839*time[i])+0.9457


    f_one_y = d*y[i-1]+e*psi[i-1]+h*delta
    f_two_y = d*(y[i-1]+grid/2)+e*(psi[i-1]+(grid/2)*f_one_y)+h*delta
    f_three_y = d*(y[i-1]+grid/2)+e*(psi[i-1]+(grid/2)*f_two_y)+h*delta
    f_four_y = d*(y[i-1]+grid)+e*(psi[i-1]+grid*f_three_y) + h*delta
    y_rk[i] = y_rk[i-1] + (grid/6)*(f_one_y+2*f_two_y+2*f_three_y+f_four_y)

    
    f_one_psi = f*y[i-1]+g*psi[i-1]+p*delta
    f_two_psi = f*(y[i-1]+grid/2)+g*(psi[i-1]+(grid/2)*f_one_psi)+p*delta
    f_three_psi = f*(y[i-1]+grid/2)+g*(psi[i-1]+(grid/2)*f_two_psi)+p*delta
    f_four_psi = f*(y[i-1]+grid)+g*(psi[i-1]+grid*f_three_psi) + p*delta
    psi_rk[i] = psi_rk[i-1] + (grid/6)*(f_one_psi+2*f_two_psi+2*f_three_psi+f_four_psi)

    lateral_acceleration[i] = d*y[i]+e*psi[i]+h*delta




'''
    y[i+1].append( y[i] + (grid_spacing)*(d*y[i]+e*psi[i]+h*delta))
    psi[i+1].append( psi(i) + (grid_spacing)*(f*y[i]+g*psi[i]+f*delta))
    time[i+1].append(time(i)+grid_spacing)
    '''

fig = plt.figure(figsize=(10,8))
ax = fig.add_subplot(321)
ax.set_title('Lateral Acceleration vs Time')
ax.set_xlabel("Time (s)")
ax.set_ylabel("Lateral Acceleration (m/s^2)")
#ax.plot(time, y, label = "Euler's Method")
#ax.plot(time, y_rk, label = "RK4")
ax.plot(time, lateral_acceleration, label = "200 km/h")
ax.legend()
plt.legend(loc="upper right")

ax = fig.add_subplot(323)
ax.set_title('Yaw Rate vs Time')
ax.set_xlabel("Time (s)")
ax.set_ylabel("Angular speed (rad/s)")
ax.plot(time, psi, label = "Euler's Method")
ax.plot(time, psi_rk, label = "RK4")
ax.plot(time, real_solution_psi, label = "Analytical Solution")
#axs[1].set_ylim(-1,3)
ax.legend()
plt.legend(loc="lower right")


ax = fig.add_subplot(322, projection='3d')
ax.set_title('3D Graph')
ax.set_xlabel("lateral speed (m/s)")
ax.set_ylabel("Angular speed (rad/s)")
ax.set_zlabel("time (s)")
surf = ax.scatter3D(lateral_acceleration,psi,time)


plt.tight_layout()
plt.show()

# fig, axs = plt.subplots(3, figsize=(6, 9), sharey=True)
# axs[0].set_title('Lateral Speed vs Time')
# axs[0].set_xlabel("Time (s)")
# axs[0].set_ylabel("Lateral Speed (m/s)")
# axs[0].plot(time, y, label = "Euler's Method")
# axs[0].plot(time, y_rk, label = "RK4")
# axs[0].plot(time, real_solution_y, label = "Analytical Solution")
# axs[0].legend()
# plt.legend(loc="upper right")


# axs[1].set_title('Yaw Rate vs Time')
# axs[1].set_xlabel("Time (s)")
# axs[1].set_ylabel("Angular speed (rad/s)")
# axs[1].plot(time, psi, label = "Euler's Method")
# axs[1].plot(time, psi_rk, label = "RK4")
# axs[1].plot(time, real_solution_psi, label = "Analytical Solution")
# #axs[1].set_ylim(-1,3)
# axs[1].legend()
# plt.legend(loc="upper right")


# axs[2] = fig.add_subplot(projection='3d')
# axs[2].set_title("3D Curve")
# axs[2].plot_surface(y,psi,time)
# axs[2].set_aspect('equal')





# plt.tight_layout()
# plt.show()
