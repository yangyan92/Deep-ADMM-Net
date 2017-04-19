function [O ,O2, DzDw1 , DzDw2 ] = xmid(I1, I2, Y, gamma, Rho, DzDy)
%% Reconstruction layer in the middle of the network.
%% This layer has two parameters: H_{l} = sum_{m=1}^{s} gamma_{l,m}B_{m};  \rho_{l}
%% This layer has two inputs: I1 = z ; I2 = \beta
%% Copyright (c) 2017 Yan Yang
%% All rights reserved.

%% network setting
config;
fN = nnconfig.FilterNumber;
fS = nnconfig.FilterSize;
pad = nnconfig.Padding;
gp = nnconfig.EnableGPU;
s = fS*fS-1;

B=filter_base( );
for i=1:fN
    H(:,:,i) = reshape( B * gamma( : ,  i) , fS ,fS);
end
HT = rot90(H,2);
load('./mask.mat')
mask = logical( ifftshift(mask) );
Denom1 = zeros(256 , 256) ; Denom1(mask) = 1 ;
[m ,n] = size(Y);
Denom2=zeros(m,m);

for k=1:fN
    prd = sqrt(Rho(k));
    Denom2 = Denom2 + abs( psf2otf (prd * H(:,:,k),[m,n])).^2;
end
Denom = Denom1+Denom2;
Denom(find(Denom == 0)) = 1e-6;
Q2=1./Denom;


if nargin == 5
    Pr=zeros(m,n);
    for i = 1:fN
        tp1 = I1(:,:,i) - I2(:,:,i);
        tp2 = imfilter(double(tp1),double(HT(:,:,i)),'same','circular','conv');
        Pr= Pr+Rho(i)*tp2;
    end
    O = real( ifft2(( fft2 ( Pr ) + Y ) .* Q2));
end

if nargin == 6
    
    if gp
        Hp = gpuArray(H);
        [Bj1p,PSp] = FBASE(B);
        Q2=gpuArray(Q2);
        DzDy=gpuArray(DzDy);
        H1p=gpuArray(HT);
        Y=gpuArray(Y);
        Ds = gpuArray(zeros(m,m));
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  O1 O2
        A = (-1)*( Q2 .* Q2);
        Trans1 =  ifft2(Q2.*fft2(DzDy));
        Trans1 = real(Trans1);
        Trans=gather(Trans1);
        PadT = padImage_g(Trans,[pad,pad],'circular');
        PadT = gpuArray(PadT);
        for i = 1:fN
            S = conv2(PadT,Hp(:,:,i),'valid');
            S = gather(S);
            O(:,:,i) = Rho(i)*S;
            O2(:,:,i) = (-1)*Rho(i)*S;
            tp5 = I1(:,:,i)-I2(:,:,i);
            aaa = padImage_g( tp5,[pad,pad],'circular');
            aaa = gpuArray(aaa);
            tp52 = conv2(aaa,H1p(:,:,i),'valid');
            Ds= Ds +Rho(i) *tp52;
        end
        Njj = fft2(Ds)+Y;
        for i = 1:fN
            PS1 = psf2otf(Hp(:,:,i), [m,m]);
            tp5 = I1(:,:,i)-I2(:,:,i);
            Padtp = padImage_g(tp5,[pad,pad],'circular');
            Padtp = gpuArray(Padtp);
            tp54 = conv2(Padtp,H1p(:,:,i),'valid');
            tp54 = fft2(tp54);
            Nii2 = abs(PS1).^2 ;
            tp1 = A.* Nii2 .*Njj;
            tp2 = Q2.*tp54;
            tp1 = ifft2(tp1);
            tp2 = ifft2(tp2);
            temp = DzDy.*(tp1+tp2);
            temp = real(gather(temp));
            DzDw2(i) = sum(temp(:));
            for  j = 1 : s
                PS3 = 2*(PSp(:,:,j).*PS1);
                tp6 = A.*PS3.*Njj;
                PS3 = ifft2(tp6);
                DS = conv2(Padtp,Bj1p(:,:,j),'valid');
                DS = fft2(DS);
                DS = Q2.*DS;
                DS = ifft2(DS);
                PfN = DzDy.*(PS3+DS);
                SSS = real(PfN*Rho(i));
                DzDw1(j,i) = sum(SSS(:));
            end
            
        end
        DzDw1=gather(DzDw1);
    else
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  O1, O2
        for k = 1:fN
            Trans =  ifft2(Q2.*fft2(DzDy));
            S = imfilter(double(Trans),double(H(:,:,k)),'same','circular','conv');
            O(:,:,k) = Rho(k)*S;
            O2(:,:,k) = (-1)*Rho(k)*S;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% DzDw1
        A=(-1)* Q2 .* Q2;
        Ds=zeros(m,m);
        for i=1:fN
            tp1 = I1(:,:,i)-I2(:,:,i);
            tp3=Rho(i) * imfilter(double(tp1) ,double(HT(:,:,i)),'same','circular','conv');
            Ds= Ds +tp3;
        end
        Njj=fft2(Ds)+Y;
        for i=1:fN
            for j=1:s
                Bj=reshape(B(:,j),fS,fS);
                Bj1=rot90(Bj,2);
                PS=psf2otf(Bj1, [m,m]);
                PS1=psf2otf(H(:,:,i), [m,m]);
                PS3=2*PS.*PS1;
                PS3=DzDy.*ifft2(Rho(i)*A.*PS3.*Njj);
                tp1=I1(:,:,i)-I2(:,:,i);
                DS= imfilter(double(tp1) ,double(Bj1),'same','circular','conv');
                DS=DzDy.*ifft2(Rho(i)*Q2.*fft2(DS));
                DzDw1(j,i)=sum(sum(PS3+DS));
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% DzDw2
        for k=1:fN
            Nii =abs(psf2otf(H(:,:,k),[m,n])).^2 ;
            tp1=I1(:,:,k) - I2(:,:,k);
            tp4 = imfilter(double(tp1),double(HT(:,:,k)),'same','circular','conv');
            tp4 =fft2(tp4);
            temp1=DzDy.*ifft2(A.*Nii.*Njj );
            temp2=DzDy.*ifft2( Q2.*tp4);
            temp=sum(temp1(:))+sum(temp2(:));
            DzDw2(k)=temp;
        end
    end
end
end
