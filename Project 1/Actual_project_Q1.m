%Euler's method
grid_spacing = 0.01;
iteration  = 1000;

initial_y = 0;
initial_psi = 0;

m = 1400;
a = 1.14;
b = 1.33;
front_stiffness = 25000;
rear_stiffness = 21000;
yaw_inertia = 2420;
velocity_x = 75;
delta = 0.1;

d = ((-(front_stiffness + rear_stiffness))/(m*velocity_x));
e = ((-(a*front_stiffness-b*rear_stiffness))/(m*velocity_x))-velocity_x;
f = ((-(a*front_stiffness-b*rear_stiffness))/(yaw_inertia*velocity_x));
g = ((-(a.^2*front_stiffness+b.^2*rear_stiffness)/(yaw_inertia*velocity_x)));


%A_Matrix = [d e;f g];
B_Matrix = [(front_stiffness/m);((a*front_stiffness)/yaw_inertia)];


y = zeros(1001,1);
psi = zeros(1001,1);
time = zeros(1001,1);

time(1) = 0;
y(1)= 0;
psi(1) = 0;

for i = 1:iteration
    y(i+1) = y(i) + (grid_spacing)*(d*y(i)+e*psi(i)+B_Matrix(1,1)*delta);
    psi(i+1) = psi(i) + (grid_spacing)*(f*y(i)+g*psi(i)+B_Matrix(2,1)*delta);
    time(i+1) = time(i)+grid_spacing;
    
end

plot(time,y);
hold on
t = linspace(0, 50, 500);
%actual_plot = -13.0964*exp(-1.9745*t) + 24.4684*exp(-0.9839*t) - 11.3720;
%plot(actual_plot);

xlabel('Time (s)');
ylabel('Velocity (m/s)');
legend ('Approximation', 'Actual');

