% Unsteady Explicit Case
clear all
close all
clc

n = 10;
nx = linspace(0,1,n);
dx = nx(2)-nx(1);
ny = nx;                                                                   % Assuming a square Grid
dy = dx;
[X,Y] = meshgrid(nx, ny);
T = 303*ones(length(nx));
T(:,1) = 400;
T(:,end) = 800;
T(1,:) = 600;
T(end,:) = 900;
cfl_alpha = -0.2;
alpha = 0.01;
dt = (cfl_alpha*(dx^2))/alpha;
k = (alpha*dt)/dx^2;
T_i = T;
t = 0;
t_max = 5;
iter = 0;
tic
while(t<t_max)
    for j=2:length(ny)-1
        for i=2:length(nx)-1
            T(i,j) =(k*(T_i(i+1,j)+T_i(i-1,j)+T_i(i,j+1)+T_i(i,j-1)-4*T_i(i,j)))+T_i(i,j);
        end 
    end
    T_i = T;
    t = t + dt;
    iter = iter + 1;
end
sim_time = toc
figure(1)
surface(nx,ny,T)
view(45,45)
%contourf(nx,ny,T)  
colorbar
saveas(figure(1),'2D_explicit.png')

