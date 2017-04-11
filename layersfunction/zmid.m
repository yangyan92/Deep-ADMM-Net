function [ O ] = zmid( p ,I1, I2, q )
% Nonlinear transform layer
% z_{l} = S_PLF(c_{l}+\beta_{l})
% The parameters are q related to the the predefined positions p;
% network setting
config;
gp = nnconfig.EnableGPU;

if nargin == 4
    I = I1+I2;
    temp = double(I);
    q = double(q);
    if gp
        temp = gpuArray(temp);
        q1 = gpuArray(q);
         p1 = gpuArray(p);
        temp2 = nnlinecu_double(p1, q1, temp);
        O = gather(temp2);
    else
        O = nnlinemex( p, q , temp);
        
    end
    
    
    
    
end

end

