function a = padImage_g(a,padSize,padding)

if ischar(padding)
    method = padding;
    padVal = [];
else
    method = 'constant';
    padVal = padding;
end

a = padarray_algo22(a, padSize, method, padVal, 'both');
