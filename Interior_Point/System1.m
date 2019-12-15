function [Ac, Bc]=System1()

    n = 12;
    m = 3;
    Ac = zeros(n,n);
    Bc = zeros(n,m);
    k = 1;

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
end