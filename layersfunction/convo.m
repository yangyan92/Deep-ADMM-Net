function [O] = convo ( I , Eta )
% Convolution layer
% This layer has a parameter: D_{l} = sum_{m=1}^{s} Eta_{l,m}B_{m};

%% network setting
config;
fN = nnconfig.FilterNumber;
fS = nnconfig.FilterSize; 
%pad = nnconfig.Padding;
%s = fS * fS - 1;

D = zeros(fS, fS, fN);
B=filter_base( ); 
for i=1:fN
    D(:,:,i) = reshape(B*Eta(:,i),fS,fS);    
end
%DT=rot90(D,2);

%%
if nargin == 2    
for i=1:fN
  O(:,:,i) = imfilter( double(I) ,double( D(:,:,i) ),'same','circular','conv');
end
end


end


    
    
    

