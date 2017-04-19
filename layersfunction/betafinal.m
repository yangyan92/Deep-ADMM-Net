function [ O, O2, O3 ,DzDw] = betafinal( I1, I2, I3, Eta, DzDy )
% The final Multiplier update layer
% O_{l} = \beta_{l}+ Eta_{l}*(c_{l}-z_{l});
%% Copyright (c) 2017 Yan Yang
%% All rights reserved.

%% network setting
config;
fN = nnconfig.FilterNumber ;

if nargin == 4
    for i = 1:1:fN
        O(:,:,i) = I1(:,:,i) +Eta(i)*( I2(:,:,i) - I3(:,:,i));
    end
end

if nargin ==5
    for i=1:fN
        O(:,:,i) = DzDy(:,:,i);
        O2(:,:,i) = Eta(i)*DzDy(:,:,i);
        O3(:,:,i) = (-1)*Eta(i)*DzDy(:,:,i);
        temp = DzDy(:,:,i).*(I2(:,:,i) - I3(:,:,i));
        DzDw(i) = sum(temp(:));
    end
end
end

