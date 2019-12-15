function [u, z0_new, mu0_new]=IP_u(Sys, T, kappa, nIter, x0, z0, mu0)

n = size(Sys.A,1);
m = size(Sys.B,2);
l = size(Sys.f,1);

% z0 = zeros(T*m+(T-1)*n, 1);
z0 = [z0(1+m+n:end,1); zeros(m+n,1)];
z = z0;
% mu0 = zeros(T*n,1);
mu = mu0;

for i =1:nIter
    g = Sys.g;
    h = Sys.h;
    b = Sys.b;
    
    g(1:m) = 2*Sys.S'*x0;
    h(1:l) = Sys.f-Sys.F1*x0;
    b(1:n,1) = Sys.A*x0;
    
    d = 1./(h-Sys.P*z);
    rd = 2*Sys.H*z+g+kappa*Sys.P'*d+Sys.C'*mu;
    rp = Sys.C*z-b;
    x = [z; mu];
    grad_A = [  2*Sys.H + kappa*Sys.P'*diag(d.^2)*Sys.P  Sys.C';
                Sys.C                           zeros(T*n,T*n)];
    grad_x = grad_A\(-[rd; rp]);
    x_new = x+1*grad_x; % This could return some meaningless values.
    z = x_new(1:T*m+(T-1)*n,1);
    mu = x_new(1+T*m+(T-1)*n:end,1);
    temp_step = 1;
    while sum(Sys.P*z>h,1)>0 
        temp_step = temp_step/2;
        x_new = x+temp_step*grad_x;
        z = x_new(1:T*m+(T-1)*n,1);
        mu = x_new(1+T*m+(T-1)*n:end,1);
    end
%     disp(sum(Sys.P*z>h,1))
%     disp(kappa*sum(-log(h-Sys.P*z),1));
end
u = z(1:m,1);
z0_new = z;
mu0_new=mu;
end