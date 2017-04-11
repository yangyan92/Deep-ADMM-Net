function [ O ] = betamid( I1, I2, I3, gamma  )
% The middle Multiplier update layer
% O_{l} = \beta_{l}+ gamma_{l}*(c_{l}-z_{l});
%% network setting
config;
fN = nnconfig.FilterNumber ;
%%
O = zeros(size(I1));
if nargin == 4
    for i = 1:fN
        O(:,:,i) = I1(:,:,i) + gamma(i) * ( I2(:,:,i) - I3(:,:,i));
    end
end
end

