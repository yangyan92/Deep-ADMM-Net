function [O, DzDw  ] = convo ( I , gamma , DzDy1 ,  DzDy2)
%% Convolution layer
%% This layer has a parameter: D_{l} = sum_{m=1}^{s} gamma_{l,m}B_{m};
%% Copyright (c) 2017 Yan Yang
%% All rights reserved.

%% network setting
config;
fN = nnconfig.FilterNumber;
fS = nnconfig.FilterSize;
pad = nnconfig.Padding;
gp = nnconfig.EnableGPU;
s=fS*fS-1;
[m,n] = size(I);

D = zeros(fS, fS, fN);
B = filter_base( );
for i = 1:fN
    D(:,:,i) = reshape(B*gamma(:,i),fS,fS);
end
DT = rot90(D,2);

if nargin == 2
    for i=1:fN
        O(:,:,i) = imfilter( double(I) ,double(D(:,:,i)),'same','circular','conv');
    end
end

if nargin == 4
    DzDy = DzDy1 + DzDy2 ;
    if gp
        Dp = gpuArray(D);
        E = gpuArray(zeros(m,m));
        D1p = gpuArray(DT);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  O
        for i = 1:fN
            d1 = padImage_g(DzDy(:,:,i),[pad,pad],'circular');
            d1 = gpuArray(d1);
            dd = conv2(d1,D1p(:,:,i),'valid');
            E = E + dd;
        end
        O = gather(E);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% DzDw
        DzDy = gpuArray(DzDy);
        Bj1p = FBASE(B,1);
        Bjp = rot90(Bj1p,2);
        DzDw = gpuArray(zeros(s,fN));
        d1 = padImage_g(I,[pad,pad],'circular');
        d1 = gpuArray(d1);
        
        
        for j = 1:fN
            for k=1:s  % m
                dd = conv2(d1,Bjp(:,:,k),'valid');
                tp = DzDy(:,:,j).* dd;
                DzDw (k,j)=sum(tp(:));
            end
        end
        
        DzDw=gather(DzDw);
        
    else
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  O
        O = zeros(m,n);
        DT = rot90(D,2);
        for i=1:fN
            O = O + imfilter(double(DzDy(:,:,i)),double(DT(:,:,i)),'same','circular','conv');
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% DzDw
        for j = 1:fN
            for k=1:s
                Bj = reshape(B(:,k),fS,fS);
                tp = DzDy(:,:,j).* imfilter(double(I),double(Bj),'same','circular','conv');
                DzDw (k,j)=sum(tp(:));
            end
        end
    end
end
end

