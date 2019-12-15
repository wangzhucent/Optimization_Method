function Sys = System_init(Ac,Bc,T,dT)
n = size(Ac,1);
m = size(Bc,2);

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
Sys.T = T;

end

