function [ O ] = betaorg( I1, I2, gamma )
% The first Multiplier updata
% O_{l} = gamma_{l}*(c_{l}-z_{l});
% network setting
config;
fN = nnconfig.FilterNumber ;
O = zeros(size(I1));
if nargin == 3
    for i = 1:fN
        O(:,:,i) = gamma(i)*(I1(:,:,i) - I2(:,:,i)) ;
    end
end

end

