function [Y, DzDw] = nnsoft(X, T, DzDy)     
% soft thresholding function 

config;
 d = nnconfig.FilterNumber;
 TH = abs(T);
 B = X >= TH;
 S = X <= -TH;

 if nargin == 2
   Y = (X - TH) .* B + (X + TH) .* S;
     DzDw = [];

else if nargin == 3
    Y = DzDy .* (B + S);
    if T >= 0
        DzDw = DzDy .* (-B + S);
    else
        DzDw = DzDy .* (B - S);
    end;
 
    else 
    error('Input arguments number not proper.');
end;
 end;

