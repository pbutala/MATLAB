function B = dec2binMat(varargin)
% dec2binMat(D) returns the binary representation of D as a Matrix.
%     D must be a non-negative integer smaller than 2^52.
%  
%     dec2binArray(D,N) produces a binary representation with at least
%     N bits.
%  
%     Example
%        dec2binMat(23) returns [1 0 1 1 1]
Bstr = dec2bin(varargin{:});
[X,N] = size(Bstr);
B = zeros(X,N);
for iN=1:N
    B(:,iN) = bin2dec(Bstr(:,iN));
end

end % end function