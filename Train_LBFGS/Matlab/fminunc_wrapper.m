function [f,g,h] = fminunc_wrapper(x,F,G,H)
% [f,g,h] = fminunc_wrapper( x, F, G, H )
% for use with Matlab's "fminunc"
f = F(x);
if nargin > 2 && nargout > 1
    g = G(x);
end
if nargin > 3 && nargout > 2
    h = H(x);
end
