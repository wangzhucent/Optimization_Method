% l = 3;
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
% x0 = zeros(n,1);
% x0 = [0.2; 0.3; -0.1; 0.1; 0.1; 0.1; zeros(6,1)];

x0 = [2; 0; 0; 0; 0; 0; zeros(6,1)];

for i = 0:T-1
    H(1+i*(m+n): (i+1)*m+i*n,1+i*(m+n):(i+1)*m+i*n) = R;    
end
for i = 1:T-1
    H(1+i*m+(i-1)*n:i*(m+n),1+i*m+(i-1)*n:i*(m+n)) = Q;
    H(1+i*m+(i-1)*n: i*(m+n), 1+i*(m+n):(i+1)*m+i*n) = S;
    H(1+i*(m+n):(i+1)*m+i*n, 1+i*m+(i-1)*n: i*(m+n)) = S';
end
g(1:m) = 2*S'*x0;

for i = 0:T-1
    P(1+i*l:(i+1)*l, 1+i*(m+n):(i+1)*m+i*n) = F2;
end
for i = 1:T-1
    P(1+i*l:(i+1)*l, 1+i*m+(i-1)*n:i*(m+n)) = F1;
end

for i = 0:T-1
    h(1+i*l:(i+1)*l) = f;
end
h(1:l) = f-F1*x0;

for i = 0:T-1
    C(1+i*n:(i+1)*n, 1+i*(m+n): (i+1)*m+i*n) = -B;
end
for i = 1:T-1
    C(1+i*n:(i+1)*n,1+i*m+(i-1)*n: i*(m+n)) = -A;
    C(1+(i-1)*n:i*n,1+i*m+(i-1)*n: i*(m+n)) = eye(n);
end

b(1:n,1) = A*x0;


nIter=10;

% z0 = zeros(T*m+(T-1)*n, 1); % Need to update for warming up 
z0 = zeros(T*m+(T-1)*n, 1);
z = z0;
mu0 = zeros(T*n,1);
mu = mu0;
z_plot =  zeros(size(z0,1), nIter);

for i =1:nIter
    d = 1./(h-P*z);
    rd = 2*H*z+g+kappa*P'*d;
    rp = C*z-b;
    x = [z; mu];
    grad_A = [  2*H + kappa*P'*diag(d.^2)*P  C';
                C                           zeros(T*n,T*n)];
    grad_x = grad_A\(-[rd; rp]);
    x = x+grad_x;
    z = x(1:T*m+(T-1)*n,1);
    mu = x(1+T*m+(T-1)*n:end,1);
    z_plot(:,i) = z;
end
for i = 1:3
    plot(1:nIter, z_plot(i,:))
    hold on
end
function y=phi(x,h,P)
    y = -ones(size(h,1))*log(h-P*x);
end
