function [f,g] = driver1(x)
% An example function, from "driver1.f" in the fortran releaes

n = length(x);
g = zeros(n,1);

f   = .25*(x(1)-1.0)^2;
for i = 2:n
    f = f + ( x(i) - (x(i-1))^2 )^2;
end
f = 4*f;

t1      = x(2) - (x(1)^2);
g(1)    = 2.0*( x(1) - 1.0 ) - 16*x(1)*t1;
for i = 2:(n-1)
    t2 = t1;
    t1 = x(i+1) - (x(i))^2;
    g(i) = 8*t2 - 16*x(i)*t1;
end
g(n)  = 8*t1;