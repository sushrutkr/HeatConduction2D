clear all
close all
clc
No_x = 50;
nx = linspace(0,1,No_x);
dx = nx(2)-nx(1);
No_y = No_x
ny = nx;                                                                   % Assuming a square Grid
dy = dx;
[x,y] = meshgrid(nx, ny);
T = 303*ones(length(nx));
T(:,1) = 400;
T(:,end) = 800;
T(1,:) = 600;
T(end,:) = 900;
k1 = 2*(dx^2+dy^2)/(dy^2);
k2 = 2*(dx^2+dy^2)/(dx^2);
T_i = T;
tolerance = 0.0001;
m_error = 9e9;
iterations = 1;
w = 1.2
method = 's';                                                              % j-jacobian g-gauss seidel s-SOR

tic
while (m_error > tolerance)

   
    for j=2:No_y-1
        for i=2:No_x-1
            if( method == 'j')                                             % Jacobi Iterative Method
                
                T(i,j) = (T_i(i+1,j)+T_i(i-1,j)+T_i(i,j+1)+T_i(i,j-1))/k1;
            
            elseif(method == 'g')                                          % Gauss Seidel Method
                
                T(i,j) = (T_i(i+1,j)+T(i-1,j)+T_i(i,j+1)+T(i,j-1))/k1;
                
            else(method == 's')                                            % Successive over relaxation
                %T(i,j) = (1/w)*((w-1)*{t_i(i,j)}+(w/k1)*(T_i(i+1,j)+T(i-1,j)+T_i(i,j+1)+T(i,j-1));
                T(i,j) = T_i(i,j)+w*(T(i-1,j)+T(i,j-1)-k1*T_i(i,j)+T_i(i+1,j)+T_i(i,j+1))/k1;
            end
        end
    end
    m_error = max(max(abs(T-T_i)));
    T_i = T;
    iterations = iterations + 1;
end
sim_time = toc
figure(1)
surface(nx,ny,T)
%contourf(nx,ny,T)
view(45,45)
