function [T, X] = Dyn_x(Sys, u0, x0, dT)
    tspan = linspace(0,dT, 10);
    odefun = @(t, x) AB(Sys, x, u0);
    [T,X] = ode45(odefun, tspan, x0);
end

function xdot = AB(Sys, x, u)
    xdot = Sys.Ac*x+Sys.Bc*u;
end