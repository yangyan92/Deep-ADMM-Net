function [ O, DzDw ] = zorg( p, I, q, DzDy1, DzDy2 )
%% The first Nonlinear transform layer
%% z_{l} = S_PLF(c_{l})
%% The parameters are q related to the the predefined positions p;
%% Copyright (c) 2017 Yan Yang
%% All rights reserved.

%% network setting
config;
gp = nnconfig.EnableGPU;
if nargin == 3
    temp = double(I);
    q = double(q);
    if gp
        temp = gpuArray(temp);
        q1 = gpuArray(q);
        p1 = gpuArray(p);
        temp2 = nnlinecu_double(p1, q1, temp);
        O = gather(temp2);
    else
        O = nnlinemex(p, q , temp);
    end
end
if nargin ==5
    DzDy = DzDy1 + DzDy2;
    xvar = double(I);
    yvar = double(DzDy);
    q = double(q);
    if gp
        xvar = gpuArray(xvar);
        yvar = gpuArray(yvar);
        q1 = gpuArray(q);
        p1 = gpuArray(p);
        [xgra, ygra] = nnlinecu_double(p1, q1, xvar, yvar);
        O = gather(xgra);
        DzDw = gather(ygra);
    else
        [O, DzDw] = nnlinemex(p, q, xvar, yvar);
    end
end
end

