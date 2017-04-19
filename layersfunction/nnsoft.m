function Y = nnsoft(X, T)   
%% Soft threshold function
%% Copyright (c) 2017 Yan Yang
%% All rights reserved.
TH = abs(T);
 B = X >= TH;
 S = X <= -TH;
 Y = (X - TH) .* B + (X + TH) .* S; 
end

