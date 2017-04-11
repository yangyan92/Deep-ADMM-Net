function [ O ] = xorg( Y, Eta, Rho )
% The first Reconstruction layer
% This layer has two parameters: H_{l} = sum_{m=1}^{s} Eta_{l,m}B_{m};  \rho_{l}


%% network setting
config;
fN = nnconfig.FilterNumber;
fS = nnconfig.FilterSize;
%pad = nnconfig.Padding;
%s = fS * fS - 1;

B=filter_base( );
H = zeros(fS, fS, fN);
for i = 1 : fN
    H(:,:,i) = reshape(B*Eta(:,i),fS,fS);
end
%DT=rot90(D,2);

load('./mask/mask_20.mat')
mask = logical( mask );
Denom1 = zeros(256 , 256) ; Denom1(mask) = 1 ;
[m ,n] = size(Y);
if nargin == 3
    Denom2 = zeros(m,m);
    for k = 1:fN
        prd = sqrt(Rho(k));
        Denom2 = Denom2 + abs( psf2otf( prd * H(:,:,k) , [m,n] ) ) .^2 ;
    end
    Denom = Denom1 + Denom2;  % diagonal matrix
    Denom(find( Denom  == 0)) = 1e-6;
    Q1 = 1./ Denom;
    O = real(ifft2(Y .* Q1)) ;
end
end











