n = 12;
m = 3;
Ac = zeros(n,n);
Bc = zeros(n,m);
dT = 0.5;

k = 1;
kappa = 0.1;

for i=1:6
    Ac(i,i+6) = 1;
end

for i = 1:6
    Ac(i+6,i) = -2*k;
end

for i = 1:5
    Ac(i+7, i) = k;
    Ac(i+6,i+1) = k;
end

Bc(7,1) = 1;
Bc(8,1) = -1;
Bc(9,2) = 1;
Bc(10,3) = 1;
Bc(11,2) = -1;
Bc(12,3) = -1;

sysc = ss(Ac,Bc,eye(n),zeros(n,m));
sysd = c2d(sysc, dT);
A = sysd.A;
B = sysd.B;

%% ------------------------------------
l=18;
F1 = zeros(l,n);
F2 = zeros(l,m);
f = zeros(l,1);

for i = 1:6
    F1(i,i) = 1;
    F1(i+6,i) = -1;
    f(i,1) = 4;
    f(i+6,1) = 4;
end

for i = 1:3
    F2(12+i,i) = 1;
    F2(15+i,i) = -1;
    f(12+i,1) = 0.5;
    f(15+i,1) = 0.5;
end

Q = eye(n);
R = eye(m);
S = zeros(n,m);

T = 30;

H = zeros(T*m+(T-1)*n, T*m+(T-1)*n);
g = zeros(T*m+(T-1)*n,1);
P = zeros(l*T, T*m+(T-1)*n);
h = zeros(l*T,1);
C = zeros(T*n, T*m+(T-1)*n);
b = zeros(T*n,1);



%% --------------------------------------
for i = 0:T-1
    H(1+i*(m+n): (i+1)*m+i*n,1+i*(m+n):(i+1)*m+i*n) = R;    
end
for i = 1:T-1
    H(1+i*m+(i-1)*n:i*(m+n),1+i*m+(i-1)*n:i*(m+n)) = Q;
    H(1+i*m+(i-1)*n: i*(m+n), 1+i*(m+n):(i+1)*m+i*n) = S;
    H(1+i*(m+n):(i+1)*m+i*n, 1+i*m+(i-1)*n: i*(m+n)) = S';
end

for i = 0:T-1
    P(1+i*l:(i+1)*l, 1+i*(m+n):(i+1)*m+i*n) = F2;
end
for i = 1:T-1
    P(1+i*l:(i+1)*l, 1+i*m+(i-1)*n:i*(m+n)) = F1;
end

for i = 0:T-1
    h(1+i*l:(i+1)*l) = f;
end

for i = 0:T-1
    C(1+i*n:(i+1)*n, 1+i*(m+n): (i+1)*m+i*n) = -B;
end
for i = 1:T-1
    C(1+i*n:(i+1)*n,1+i*m+(i-1)*n: i*(m+n)) = -A;
    C(1+(i-1)*n:i*n,1+i*m+(i-1)*n: i*(m+n)) = eye(n);
end

%% Initialization of System
Sys.Ac = Ac;
Sys.Bc = Bc;
Sys.A = A;
Sys.B = B;
Sys.F1 = F1;
Sys.F2 = F2;
Sys.f = f;
Sys.Q = Q;
Sys.S = S;
Sys.H = H;
Sys.h = h;
Sys.g = g;
Sys.P = P;
Sys.C = C;
Sys.b = b;

Sys.x0 = x0;

z0 = zeros(T*m+(T-1)*n, 1);
mu0 = zeros(T*n,1);
x0 = [2; 0; 0; 0; 0; 0; zeros(6,1)];

u0 = [1;0; 0];

nT = 100;
x_plot = zeros(size(x0,1),nT);
for i = 1: nT
    [T,X] = Dyn_x(Sys, u0,x0, dT);
    x0 = X(end,:)';
    x_plot(:,i) = x0;
end
figure(1)
for i = 1:2
    plot(x_plot(i,:))
    hold on
end 




