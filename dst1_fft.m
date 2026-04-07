function Y = dst1_fft(X, dim)
%DST1_FFT Discrete sine transform of type-I computed by FFT.
%   Y = DST1_FFT(X) applies DST-I along the first non-singleton dimension.
%   Y = DST1_FFT(X, DIM) applies DST-I along dimension DIM.
%
%   This implementation is self-contained and does not depend on any toolbox.
%   The inverse relation is:
%       X = 2/(n+1) * DST1_FFT(Y, DIM)
%   where n = size(X, DIM).

if nargin < 2 || isempty(dim)
    dim = find(size(X) ~= 1, 1, 'first');
    if isempty(dim)
        dim = 1;
    end
end

sz = size(X);
n = sz(dim);
perm = [dim, 1:dim-1, dim+1:ndims(X)];
Xp = permute(X, perm);
szp = size(Xp);
Xp = reshape(Xp, n, []);

Yext = zeros(2 * (n + 1), size(Xp, 2));
Yext(2:n+1, :) = Xp;
Yext(n+3:end, :) = -flipud(Xp);
Yfft = fft(Yext, [], 1);
Yp = -imag(Yfft(2:n+1, :)) / 2;

Yp = reshape(Yp, szp);
Y = ipermute(Yp, perm);
end
