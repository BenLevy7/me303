%Euler's method
%The smaller the grid spacing, the more accurate the result
syms x
grid_spacing = 0.01;
%y_next = 0;
y_current = 1; %declared based on initial conditions
iteration = 100;

for i = 1:iteration
    y_next = y_current + (grid_spacing)*(y_current);
    y_current = y_next;
    
end 

disp(y_current);