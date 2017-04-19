function Y = rnnloss(X, I, DzDy )
%% rnnloss: calculate the NMSE of restored image and original image
%% X: reconstructed image of size m*n;
%% I: ground-truth image of size m*n;

X = double(X);
I = double(I);
B=norm(I,'fro');

if nargin == 2
 S = X - I ;
 Y = norm(S,'fro') / B ;
elseif nargin ==3
 S = X - I ;
 Y1 = norm(S,'fro') ;   
 Y = S /(B*Y1);
else
    error('Input arguments number not proper.');
end;
end