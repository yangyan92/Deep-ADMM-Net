function aIdx = getPaddingIndices22(aSize,padSize,method,direction)
%getPaddingIndices is used by padarray and blockproc. 
%   Computes padding indices of input image.  This is function is used to
%   handle padding of in-memory images (via padarray) as well as
%   arbitrarily large images (via blockproc).
%
%   aSize : result of size(I) where I is the image to be padded
%   padSize : padding amount in each dimension.  
%             numel(padSize) can be greater than numel(aSize)
%   method : X or a 'string' padding method
%   direction : pre, post, or both.
%
%   See the help for padarray for additional information.

% Copyright 2010 The MathWorks, Inc.

% make sure we have enough image dims for the requested padding
if numel(padSize) > numel(aSize)
    singleton_dims = numel(padSize) - numel(aSize);
    aSize = [aSize ones(1,singleton_dims)];
end

switch method
    case 'circular'
        aIdx = CircularPad(aSize, padSize, direction);
    case 'symmetric'
        aIdx = SymmetricPad(aSize, padSize, direction);
    case 'replicate' 
        aIdx = ReplicatePad(aSize, padSize, direction);
end


%%%
%%% CircularPad
%%%
function idx = CircularPad(aSize, padSize, direction)

numDims = numel(padSize);

% Form index vectors to subsasgn input array into output array.
% Also compute the size of the output array.
idx   = cell(1,numDims);
for k = 1:numDims
    M = aSize(k);
    dimNums = uint32(1:M);
    p = padSize(k);
    
    switch direction
        case 'pre'
            idx{k}   = dimNums(mod(-p:M-1, M) + 1);
            
        case 'post'
            idx{k}   = dimNums(mod(0:M+p-1, M) + 1);
            
        case 'both'
            idx{k}   = dimNums(mod(-p:M+p-1, M) + 1);
            
    end
end


%%%
%%% SymmetricPad
%%%
function idx = SymmetricPad(aSize, padSize, direction)

numDims = numel(padSize);

% Form index vectors to subsasgn input array into output array.
% Also compute the size of the output array.
idx   = cell(1,numDims);
for k = 1:numDims
    M = aSize(k);
    dimNums = uint32([1:M M:-1:1]);
    p = padSize(k);
    
    switch direction
        case 'pre'
            idx{k}   = dimNums(mod(-p:M-1, 2*M) + 1);
            
        case 'post'
            idx{k}   = dimNums(mod(0:M+p-1, 2*M) + 1);
            
        case 'both'
            idx{k}   = dimNums(mod(-p:M+p-1, 2*M) + 1);
    end
end


%%%
%%% ReplicatePad
%%%
function idx = ReplicatePad(aSize, padSize, direction)

numDims = numel(padSize);

% Form index vectors to subsasgn input array into output array.
% Also compute the size of the output array.
idx   = cell(1,numDims);
for k = 1:numDims
    M = aSize(k);
    p = padSize(k);
    onesVector = uint32(ones(1,p));
    
    switch direction
        case 'pre'
            idx{k}   = [onesVector 1:M];
            
        case 'post'
            idx{k}   = [1:M M*onesVector];
            
        case 'both'
            idx{k}   = [onesVector 1:M M*onesVector];
    end
end

