function [u,z] = IP_quadprog(Sys, T, kappa, nIter, x0, z0, mu0)
    n = size(Sys.A,1);
    m = size(Sys.B,2);
    l = size(Sys.f,1);
    g = Sys.g;
    h = Sys.h;
    b = Sys.b;
    
    g(1:m) = 2*Sys.S'*x0;
    h(1:l) = Sys.f-Sys.F1*x0;
    b(1:n,1) = Sys.A*x0;
    
    z0 = [z0(1+m+n:end,1); zeros(m+n,1)];
    z = z0;
    
    options = optimoptions('quadprog', 'OptimalityTolerance', 0,'StepTolerance',0,'MaxIterations',20);
    z = quadprog(Sys.H*2, g, Sys.P, h, Sys.C, b, [], [], z0, options);
    u = z(1:m,1);
end