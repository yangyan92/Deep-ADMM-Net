function [ O, O2 ,DzDw ] = betaorg( I1, I2, Eta, DzDy1, DzDy2, DzDy3  )
%% The first Multiplier updata
%% O_{l} = Eta_{l}*(c_{l}-z_{l});
%% Copyright (c) 2017 Yan Yang
%% All rights reserved.

%% network setting
config;
fN = nnconfig.FilterNumber;

if nargin == 3
    for i = 1:1:fN
        O(:,:,i) = Eta(i)*(I1(:,:,i) - I2(:,:,i)) ;
    end
end

if nargin == 6
    DzDy = DzDy1 + DzDy2 +DzDy3;
    for i =1:fN
        O(:,:,i) = Eta(i) * DzDy(:,:,i);
        O2 (:,:,i) = (-1) * Eta(i) * DzDy(:,:,i);
        temp = DzDy(:,:,i).*(I1(:,:,i) - I2(:,:,i));
        DzDw(i) = sum(temp(:));
    end
end


end

