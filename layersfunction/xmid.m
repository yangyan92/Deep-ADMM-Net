function [O ] = xmid(I1, I2, Y, Eta , Rho)
% Reconstruction layer in the middle of the network.
% This layer has two parameters: H_{l} = sum_{m=1}^{s} Eta_{l,m}B_{m};  \rho_{l}
% This layer has two inputs: I1 = z ; I2 = \beta

%% network setting
config;
fN = nnconfig.FilterNumber;
fS = nnconfig.FilterSize;
%pad = nnconfig.Padding;
%s = fS * fS - 1;


B=filter_base( );
for i=1:fN
    H(:,:,i) = reshape( B * Eta( : ,  i) , fS ,fS);
end
HT = rot90(H,2);

load('./mask/mask_20.mat')
mask = logical( ifftshift(mask) );
Denom1 = zeros(256 , 256) ; Denom1(mask) = 1 ;
[m ,n] = size(Y);

%%
if nargin == 5
    Denom2=zeros(m,m);
    for k=1:fN
       prd = sqrt(Rho(k));
       Denom2 = Denom2 + abs( psf2otf (prd * H(:,:,k),[m,n])).^2;
    end
    Denom = Denom1+Denom2;
    Denom(find(Denom == 0)) = 1e-6;
    Q2=1./Denom;
    A=zeros(m,n);
    for i = 1:fN
        tp1 = I1(:,:,i) - I2(:,:,i); %Zh
        tp2 = imfilter(double(tp1),double(HT(:,:,i)),'same','circular','conv'); %Zh2
        A= A+Rho(i)*tp2;
    end
    O = real( ifft2(( fft2 ( A ) + Y ) .* Q2));
    
end


end
