function [Ac,Bc] = System2()
%{ 
    This is a system of 
    Ac =    0   0   0   1   0   0
            0   0   0   0   1   0
            0   0   0   0   0   1
            0   0   0   0   0   0
            0   0   0   0   0   0
            0   0   0   0   0   0
%}
    n = 6;
    m = 3;
    J1 = 120;
    J2 = 100;
    J3 = 80;
    Ac = zeros(n,n);
    Bc = zeros(n,m);
    for i=1:3
        Ac(i,i+3) = 1;
    end
    Bc(4,1) = 1/J1;
    Bc(5,2) = 1/J2;
    Bc(6,3) = 1/J3;
end