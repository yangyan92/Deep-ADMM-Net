function [Bj1p,PSp] = FBASE( B,a)
config;
fS = nnconfig.FilterSize;
s=fS*fS-1;
size = nnconfig.ImageSize ;

if nargin == 1
    for j = 1:s
        Bj = reshape(B(:,j),fS,fS);
        Bj1(:,:,j) = rot90(Bj,2);
        PS(:,:,j) = psf2otf(Bj1(:,:,j), size);
    end
    PSp =gpuArray(PS);
    Bj1p = gpuArray(Bj1);
end
if nargin == 2
    for j = 1:s
        Bj = reshape(B(:,j),fS,fS);
        Bj1(:,:,j) = rot90(Bj,2);
    end
    Bj1p = gpuArray(Bj1)  ;
end
end

