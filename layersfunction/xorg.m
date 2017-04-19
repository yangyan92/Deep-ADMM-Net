function [ O , DzDw1 , DzDw2 ] = xorg( Y, gamma, Rho,  DzDy )
%% The first Reconstruction layer
%% This layer has two parameters: H_{l} = sum_{m=1}^{s} gamma_{l,m}B_{m};  \rho_{l}
%% Copyright (c) 2017 Yan Yang
%% All rights reserved.

%% network setting
config;
fN = nnconfig.FilterNumber;
fS = nnconfig.FilterSize;
gp = nnconfig.EnableGPU;
s = fS*fS-1;
[m ,n] = size(Y);
%%
B = filter_base( );
H = zeros(fS, fS, fN);
for i = 1 : fN
    H(:,:,i) = reshape(B*gamma(:,i),fS,fS);
end
HT = rot90(H,2);
load('./mask.mat')
mask = logical( ifftshift(mask) );
Denom1 = zeros(m , n) ; Denom1(mask) = 1 ;
Denom2 = zeros(m,m);
for k = 1:fN
    prd = sqrt(Rho(k));
    Denom2 = Denom2 + abs( psf2otf( prd * H(:,:,k), [m,n] ) ) .^2 ;
end
Denom = Denom1 + Denom2;  % diagonal matrix
Denom(find( Denom  == 0)) = 1e-6;
Q1 = 1./ Denom;

%% The forward propagation
if nargin == 3
    O = real(ifft2(Y .* Q1)) ;
end

%% The backward propagation
if nargin == 4
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%O
    O = 1;
    if gp
        Hp=gpuArray(H);
        H1p=gpuArray(HT);
        [Bj1p,PSp]=FBASE(B);
        Bjp = rot90(Bj1p,2);
        Q1 = gpuArray(Q1);
        DzDy = gpuArray(DzDy);
        Y = gpuArray(Y);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%DzDw1  DzDw2
        A = -1*Q1.*Q1;
        for i = 1:fN
            PS1 = psf2otf(Hp(:,:,i), [m,m]);
            tp1 = abs(PS1).^2 ;
            tp2 = A.*tp1.*Y;
            tp2 = ifft2(tp2);
            tp2 = real(tp2);
            temp = DzDy.*tp2;
            DzDw2(i) = sum(temp(:));
            for j = 1:s
                PS3 = 2*(PSp(:,:,j).*PS1.*Y);
                tp3 = Rho(i)*(A.*PS3);
                tp3 = real(ifft2(tp3));
                PS4 = DzDy.*tp3;
                DzDw1(j,i) = sum(sum(PS4));
            end
        end
        DzDw2=gather(DzDw2);
        DzDw1=gather(DzDw1);               
    else        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% DzDw1
        A = (-1) * Q1 .*Q1 ;
        for i = 1:fN
            for j=1:s
                Bj = reshape(B(:,j),fS,fS);
                Bj1 = rot90(Bj,2);
                PS = psf2otf(Bj1, [m,m]);
                PS1 = psf2otf(H(:,:,i), [m,m]);
                PS3 = 2*PS.*PS1.*Y;
                PS4 = DzDy.*ifft2(Rho(i)*A.*PS3);
                DzDw1(j,i) = sum(sum(PS4));
            end
        end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% DzDw2
        for k=1:fN
            APF = abs(psf2otf(H(:,:,k),[m,n])).^2 ;
            temp = DzDy.*ifft2(A.*APF.*Y);
            DzDw2(k) = sum(temp(:));
        end
    end
end
end
