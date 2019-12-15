clear all
clc
% dT = 0.5;
dT = 2;
kappa = 0.1;
T = 50;

% [Ac, Bc] =System1();
[Ac, Bc] = System2();
n = size(Ac,1);
m = size(Bc,2);

% Sys = System_init(Ac,Bc,T,dT);
Sys = System_init(Ac,Bc,T,dT);
z0 = zeros(T*m+(T-1)*n, 1);
z = z0;
mu0 = zeros(T*n,1);
mu = mu0;
% x0 = [2; -2; 0; 0; 0; 0; zeros(6,1)];
%  x0 = [2; 0; 0; 0; 0; 0; zeros(6,1)];
x0 = [-0.4; -0.8; 1.2; -0.02; -0.02; 0.02];

x = x0;
nIter = 20;

t_plot =[];
x_plot =[];
u_plot = [];
tic

for i = 1:50
    i
    
    [u, z, mu]=IP_u_faster(Sys, T, kappa, nIter, x, z, mu);
%     [u, z] = IP_quadprog(Sys, T, kappa, nIter, x, z, mu);
    [temp_T, temp_X] = Dyn_x(Sys, u, x, dT);
    x = temp_X(end,:)';
    t_plot = [t_plot; temp_T+(i-1)*dT];
    x_plot = [x_plot; temp_X];
    u_plot = [u_plot; u'];
end
ntime = toc;
disp(ntime)
figure(1)
for i = 1:6
    plot(t_plot,x_plot(:,i))
    hold on
end

figure(2)
for i = 1:3
    stairs(u_plot(:,i))
    hold on
end





