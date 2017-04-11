function [ O ] = zorg( p , I, q )
% The first Nonlinear transform layer
% z_{l} = S_PLF(c_{l})
% The parameters are q related to the the predefined positions p;


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


end

