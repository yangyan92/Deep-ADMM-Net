function b = padarray_algo22(a, padSize, method, padVal, direction)
%PADARRAY_ALGO Pad array.
%   B = PADARRAY_AGLO(A,PADSIZE,METHOD,PADVAL,DIRECTION) internal helper
%   function for PADARRAY, which performs no input validation.  See the
%   help for PADARRAY for the description of input arguments, class
%   support, and examples.

%   Copyright 2014 The MathWorks, Inc.

if isempty(a)
    
    numDims = numel(padSize);
    sizeB = zeros(1,numDims);
    
    for k = 1: numDims
        % treat empty matrix similar for any method
        if strcmp(direction,'both')
            sizeB(k) = size(a,k) + 2*padSize(k);
        else
            sizeB(k) = size(a,k) + padSize(k);
        end
    end
    
    b = mkconstarray(class(a), padVal, sizeB);
    
elseif strcmpi(method,'constant')
    
    % constant value padding with padVal
    b = ConstantPad(a, padSize, padVal, direction);
else
    
    % compute indices then index into input image
    aSize = size(a);
    aIdx = getPaddingIndices22(aSize,padSize,method,direction);
    b = a(aIdx{:});
end

if islogical(a)
    b = logical(b);
end

%%%
%%% ConstantPad
%%%
function b = ConstantPad(a, padSize, padVal, direction)

numDims = numel(padSize);

% Form index vectors to subsasgn input array into output array.
% Also compute the size of the output array.
idx   = cell(1,numDims);
sizeB = zeros(1,numDims);
for k = 1:numDims
    M = size(a,k);
    switch direction
        case 'pre'
            idx{k}   = (1:M) + padSize(k);
            sizeB(k) = M + padSize(k);
            
        case 'post'
            idx{k}   = 1:M;
            sizeB(k) = M + padSize(k);
            
        case 'both'
            idx{k}   = (1:M) + padSize(k);
            sizeB(k) = M + 2*padSize(k);
    end
end

% Initialize output array with the padding value.  Make sure the
% output array is the same type as the input.
b         = mkconstarray(class(a), padVal, sizeB);
b(idx{:}) = a;