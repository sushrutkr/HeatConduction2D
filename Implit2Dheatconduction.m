%Implicit method
clear all
close all
clc

nx = 10;
Nx = linspace(0,1,nx);
dx = Nx(2)-Nx(1);
Ny = Nx;
ny = nx;                                                                   % Assuming a square Grid
dy = dx;
[x,y] = meshgrid(nx, ny);
T = 303*ones(nx);
T(:,1) = 400;
T(:,end) = 800;
T(1,:) = 900;
T(end,:) = 600;
alpha = 0.01;
cfl_alpha = 0.5;
dt = (cfl_alpha*(dx^2))/alpha;
k1 = alpha*(dt/(2*(dx^2)));
k2 = alpha*(dt/(2*(dy^2)));
T_i = T;
n = nx^2;
A = zeros(n);
B = zeros(n,1);
t_max = 5;
t = 0;
iter = 0;
tic

while (t_max > t)

            %Left Boundaries
            j=1;
            for i = 1:nx
                 A((j-1)*ny+i,(j-1)*ny+i) = 1;
                 B((j-1)*ny+i) = 400;
            end
            
            %Right Boundaries
            j = ny;
            for i = 1:ny
                 A((j-1)*ny+i,(j-1)*ny+i) = 1;
                 B((j-1)*ny+i) = 800;
            end
            
            %Top Boundaries
            i = 1;
            for j = 1:ny
                A((j-1)*ny+i,(j-1)*ny+i) = 1;
                B((j-1)*ny+i) = 600;
            end
            
            %Bottom Boundaries
            i=nx;
            for j=1:ny
                 A((j-1)*ny+i,(j-1)*ny+i) = 1;
                 B((j-1)*ny+i) = 900;
            end
             %Innernodes
                 for j = 2:ny-1
                     for i =2:nx-1
                        % Inner Nodes
                        A((j-1)*nx+i,(j-1)*nx+i)   = 1+4*k1;
                        A((j-1)*nx+i,(j-1)*nx+i-1) = -k1;
                        A((j-1)*nx+i,(j-1)*nx+i+1) = -k1;
                        A((j-1)*nx+i,(j-2)*nx+i)   = -k2;
                        A((j-1)*nx+i,(j)*nx+i)     = -k2;
                        B((j-1)*nx+i)              = T_i(i,j); 
                     end
                 end
                 
            T_1D = inv(A)*(B);
            T = reshape(T_1D,nx,ny);
            T_i = T;
            t = t + dt;
            iter = iter +1;
            
end  
sim_time = toc

figure(1)
surface(Nx,Ny,T)
view(45,45)            
           