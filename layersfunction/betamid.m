function [ O, O2, O3,DzDw ] = betamid( I1, I2, I3, Eta, DzDy1, DzDy2, DzDy3  )
%% The middle Multiplier update layer
%% O_{l} = \beta_{l}+ Eta_{l}*(c_{l}-z_{l});
%% Copyright (c) 2017 Yan Yang
%% All rights reserved.

%% network setting 
config;
fN = nnconfig.FilterNumber ; 

if nargin == 4  
     for i = 1:1:fN
      O(:,:,i) = I1(:,:,i) + Eta(i)*( I2(:,:,i) - I3(:,:,i));
     end
end

if nargin == 7
    for i = 1:fN
    DzDy = DzDy1 + DzDy2 +DzDy3;
    O(:,:,i) = DzDy(:,:,i);
    O2(:,:,i) = Eta(i)*DzDy(:,:,i);
    O3(:,:,i) = (-1)*Eta(i)*DzDy(:,:,i);
    temp = DzDy(:,:,i).*(I2(:,:,i) - I3(:,:,i));
    DzDw(i) = sum(temp(:));
    end
end


end

