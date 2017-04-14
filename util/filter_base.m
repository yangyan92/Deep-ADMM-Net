function B = filter_base( )
%%%%%%%%%%%%%%%%%%%%DCT base
 config;
 fS = nnconfig.FilterSize ;
 fN = nnconfig.FilterNumber;
 fS_sqrt = fS^2;

 DCT = dctmtx(fS);
 DCT = kron(DCT, DCT);
 B = zeros(fS_sqrt, fN);
 for i = 2 : fS_sqrt
 B(:, i-1) = DCT(i, :);
 end


end

