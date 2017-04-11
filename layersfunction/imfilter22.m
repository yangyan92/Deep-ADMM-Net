function b = imfilter(varargin)
%IMFILTER N-D filtering of multidimensional images.
%   B = IMFILTER(A,H) filters the multidimensional array A with the
%   multidimensional filter H.  A can be logical or it can be a 
%   nonsparse numeric array of any class and dimension.  The result, 
%   B, has the same size and class as A.
%
%   Each element of the output, B, is computed using double-precision
%   floating point.  If A is an integer or logical array, then output 
%   elements that exceed the range of the given type are truncated, 
%   and fractional values are rounded.
%
%   B = IMFILTER(A,H,OPTION1,OPTION2,...) performs multidimensional
%   filtering according to the specified options.  Option arguments can
%   have the following values:
%
%   - Boundary options
%
%       X            Input array values outside the bounds of the array
%                    are implicitly assumed to have the value X.  When no
%                    boundary option is specified, IMFILTER uses X = 0.
%
%       'symmetric'  Input array values outside the bounds of the array
%                    are computed by mirror-reflecting the array across
%                    the array border.
%
%       'replicate'  Input array values outside the bounds of the array
%                    are assumed to equal the nearest array border
%                    value.
%
%       'circular'   Input array values outside the bounds of the array
%                    are computed by implicitly assuming the input array
%                    is periodic.
%
%   - Output size options
%     (Output size options for IMFILTER are analogous to the SHAPE option
%     in the functions CONV2 and FILTER2.)
%
%       'same'       The output array is the same size as the input
%                    array.  This is the default behavior when no output
%                    size options are specified.
%
%       'full'       The output array is the full filtered result, and so
%                    is larger than the input array.
%
%   - Correlation and convolution
%
%       'corr'       IMFILTER performs multidimensional filtering using
%                    correlation, which is the same way that FILTER2
%                    performs filtering.  When no correlation or
%                    convolution option is specified, IMFILTER uses
%                    correlation.
%
%       'conv'       IMFILTER performs multidimensional filtering using
%                    convolution.
%
%   Notes
%   -----
%   This function may take advantage of hardware optimization for datatypes
%   uint8, uint16, int16, single, and double to run faster.
%
%   Example 
%   -------------
%       originalRGB = imread('peppers.png'); 
%       h = fspecial('motion',50,45); 
%       filteredRGB = imfilter(originalRGB,h); 
%       figure, imshow(originalRGB), figure, imshow(filteredRGB)
%       boundaryReplicateRGB = imfilter(originalRGB,h,'replicate'); 
%       figure, imshow(boundaryReplicateRGB)
%
%   See also FSPECIAL, CONV2, CONVN, FILTER2. 

%   Copyright 1993-2013 The MathWorks, Inc.

% Testing notes
% Syntaxes
% --------
% B = imfilter(A,H)
% B = imfilter(A,H,Option1, Option2,...)
%
% A:       numeric, full, N-D array.  May not be uint64 or int64 class. 
%          May be empty. May contain Infs and Nans. May be complex. Required.
%         
% H:       double, full, N-D array.  May be empty. May contain Infs and Nans.
%          May be complex. Required.
%
% A and H are not required to have the same number of dimensions. 
%
% OptionN  string or a scalar number. Not case sensitive. Optional.  An
%          error if not recognized.  While there may be up to three options
%          specified, this is left unchecked and the last option specified
%          is used.  Conflicting or inconsistent options are not checked.
%
%        A choice between these options for boundary options
%        'Symmetric' 
%        'Replicate'
%        'Circular'
%         Scalar #  - Default to zero.
%       A choice between these strings for output options
%        'Full'
%        'Same'  - default
%       A choice between these strings for functionality options
%        'Conv' 
%        'Corr'  - default       
%
% B:   N-D array the same class as A.  If the 'Same' output option was
%      specified, B is the same size as A.  If the 'Full' output option was
%      specified the size of B is size(A)+size(H)-1, remembering that if
%      size(A)~=size(B) then the missing dimensions have a size of 1.
%
% 
% IMFILTER should use a significantly less amount of memory than CONVN. 

% MATLAB Compiler pragma: iptgetpref is indirectly invoked by the code that
% loads the Intel IPP library.
%#function iptgetpref

[a, h, boundary, sameSize, convMode] = parse_inputs(varargin{:});

[finalSize, pad] = computeSizes(a, h, sameSize);

%Empty Inputs
% 'Same' output then size(b) = size(a)
% 'Full' output then size(b) = size(h)+size(a)-1 
if isempty(a)
  
  b = handleEmptyImage(a, sameSize, finalSize);
  return
  
elseif isempty(h)
    
  b = handleEmptyFilter(a, sameSize, finalSize);
  return
  
end

% Pad input based on dimensions of filter kernel.
a = padImage(a,pad,boundary);

% Separate real and imaginary parts of the filter (h) in MATLAB and
% filter imaginary and real parts of the image (a) in the mex code. 
if (isSeparable(a, h))
    
  % extract the components of the separable filter
  [u,s,v] = svd(h);
  s = diag(s);
  hcol = u(:,1) * sqrt(s(1));
  hrow = v(:,1)' * sqrt(s(1));
  
  % intermediate results should be stored in doubles in order to
  % maintain sufficient precision
  class_of_a = class(a);
  if ~isa(a,'double')
    change_class = true;
    a = double(a);
  else
    change_class = false;
  end
  
  % apply the first component of the separable filter (hrow)
  out_size_row = [size(a,1) finalSize(2:end)];                        
  start = [0 pad(2:end)];
  b_tmp = filterPartOrWhole(a, out_size_row, hrow, start, sameSize, convMode);
  
  % apply the other component of the separable filter (hcol)
  start = [pad(1) 0 pad(3:end)];
  b = filterPartOrWhole(b_tmp, finalSize, hcol, start, sameSize, convMode);
  
  if change_class
      b = cast(b, class_of_a);
  end
  
else % non-separable filter case
    a1=gpuArray(a);
    finalSize1=gpuArray(finalSize);
    h1=gpuArray(h);
    pad1=gpuArray(pad);
    sameSize1=gpuArray(sameSize);
  convMode1=gpuArray(convMode);
  b1 = filterPartOrWhole(a1, finalSize1, h1, pad1, sameSize1, convMode1);
  b=gather(b1);
    
end

%======================================================================

%--------------------------------------------------------------
function [a, h, boundary, sameSize, convMode] = parse_inputs(a, h, varargin)

narginchk(2,5);

validateattributes(a,{'numeric' 'logical'},{'nonsparse'},mfilename,'A',1);
validateattributes(h,{'double'},{'nonsparse'},mfilename,'H',2);

%Assign defaults
boundary = 0;  %Scalar value of zero
output = 'same';
do_fcn = 'corr';

allStrings = {'replicate', 'symmetric', 'circular', 'conv', 'corr', ...
              'full','same'};

for k = 1:length(varargin)
  if ischar(varargin{k})
    string = validatestring(varargin{k}, allStrings,...
                          mfilename, 'OPTION',k+2);
    switch string
     case {'replicate', 'symmetric', 'circular'}
      boundary = string;
     case {'full','same'}
      output = string;
     case {'conv','corr'}
      do_fcn = string;
    end
  else
    validateattributes(varargin{k},{'numeric'},{'nonsparse'},mfilename,'OPTION',k+2);
    boundary = varargin{k};
  end %else
end

sameSize = strcmp(output,'same');

convMode = strcmp(do_fcn,'conv');

%--------------------------------------------------------------
function separable = isSeparable(a, h)

% check for filter separability 
sep_threshold = getSeparableFilterThreshold(class(a));

if ((numel(h) >= sep_threshold) && ...
    (ismatrix(h)) && ...
    all(size(h) ~= 1) && ...
    all(isfinite(h(:))))

  [~,s,~] = svd(h);
  s = diag(s);
  tol = length(h) * max(s) * eps;
  rank = sum(s > tol);
  
  if (rank == 1)
    separable = true;
  else
    separable = false;
  end
  
else
    
  separable = false;
    
end

%--------------------------------------------------------------
function b = handleEmptyImage(a, sameSize, im_size)

if (sameSize)

  b = a;
    
else
    
  if all(im_size >= 0)
      
    b = zeros(im_size, class(a));
      
  else
      
    error(message('images:imfilter:negativeDimensionBadSizeB'))
      
  end
    
end

%--------------------------------------------------------------
function b = handleEmptyFilter(a, sameSize, im_size)

if (sameSize)
  
    b = zeros(size(a), class(a));
  
else
  
  if all(im_size>=0)
    
    b = zeros(im_size, class(a));
    
  else
    
    error(message('images:imfilter:negativeDimensionBadSizeB'))
    
  end
  
end

%--------------------------------------------------------------
function ippFlag = useIPPL(a,outSize,h,nonzero_h)

% Querying the Java preference for IPP use is computationally expensive,
% particularly in cases where imfilter is called in a loop. As a
% performance optimization, store the preference as persistent state. This
% state is reset whenever iptsetpref('UseIPPL',___) is called.
persistent prefFlag;
needToQueryJavaPreference = isempty(prefFlag);
if needToQueryJavaPreference
    prefFlag = iptgetpref('UseIPPL');
end

%We are disabling the use of IPP for double precision inputs on win32.
if ~isImageIPPFilterType(class(a))
    ippFlag = false;
    return;
end

if numel(nonzero_h)/numel(h) > 0.05
    densityFlag = true;
else
    densityFlag = false;
end

hDimsFlag = ismatrix(h);

% Determine if the image is big depending on datatype
tooBig = isImageTooBigForIPPFilter(a, outSize);

ippFlag = prefFlag && densityFlag && hDimsFlag && (~tooBig);

%--------------------------------------------------------------
function [finalSize, pad] = computeSizes(a, h, sameSize)

rank_a = ndims(a);
rank_h = ndims(h);

% Pad dimensions with ones if filter and image rank are different
size_h = [size(h) ones(1,rank_a-rank_h)];
size_a = [size(a) ones(1,rank_h-rank_a)];

if (sameSize)
  %Same output
  finalSize = size_a;

  %Calculate the number of pad pixels
  filter_center = floor((size_h + 1)/2);
  pad = size_h - filter_center;
else
  %Full output
  finalSize = size_a+size_h-1;
  pad = size_h - 1;
end                        

%--------------------------------------------------------------
function result = filterPartOrWhole(a, outSize, h, start, sameSize, convMode)

% Create connectivity matrix.  Only use nonzero values of the filter.
conn = h~=0;
nonzero_h = h(conn);

ippFlag  = useIPPL(a,outSize,h,nonzero_h);

if (isreal(h))
  
  result = imfilter_mex(a, outSize, h, nonzero_h,...
                        conn, start, sameSize, convMode, ippFlag);
  
else
  
  b1 = imfilter_mex(a, outSize, real(h), real(nonzero_h),...
                    conn, start, sameSize, convMode, ippFlag);
  
  b2 = imfilter_mex(a, outSize, imag(h), imag(nonzero_h), ...
                    conn, start, sameSize, convMode, ippFlag);
  
  if isreal(a)
    
    % b1 and b2 will always be real; result will always be complex
    result = complex(b1, b2);
    
  else
    
    % b1 and/or b2 may be complex;
    % result will always be complex
    result = complex(real(b1) - imag(b2),...
                     imag(b1) + real(b2));
    
  end
  
end

%--------------------------------------------------------------
function a = padImage(a,padSize,padding)

if ischar(padding)
    method = padding;
    padVal = [];
else
    method = 'constant';
    padVal = padding;
end

a = padarray_algo(a, padSize, method, padVal, 'both');
