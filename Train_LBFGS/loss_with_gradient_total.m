function [ loss, grad ] = loss_with_gradient_total( weights)
%% computing the loss and the gradients of the all training samples.
%% Copyright (c) 2017 Yan Yang 
%% All rights reserved.
%%
net = weiTOnet(weights);
loss = 0;
config;
TN = nnconfig.TrainNumber;

if nargout == 1
    for i = 1:1:TN
        % get training data
        data = getMData(i);  
         l = loss_with_gradient_single(data, net);
        loss = loss + l;
    end;
    loss = loss / TN;
    
elseif nargout == 2
    grad_length = length(weights);
    grad = zeros(grad_length,1);
    for i = 1:1:TN
        data = getMData(i);
        [l, g] = loss_with_gradient_single(data, net);
        loss = loss + l;      
        grad = grad + g;
    end;
    loss = loss / TN;
    grad = grad / TN;
else
    error('Invalid out put number.');
end;

end

