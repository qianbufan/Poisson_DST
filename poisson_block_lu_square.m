function [x, y, U, info] = poisson_block_lu_square(N, f_handle)

if N < 2 || floor(N) ~= N
    error('N must be an integer greater than or equal to 2.');
end

n = N - 1;
h = 1 / N;
x = linspace(0, 1, N + 1).';
y = linspace(0, 1, N + 1);
t = h * (1:n);
[X, Y] = ndgrid(t, t);
F = h^2 * double(f_handle(X, Y));

D = 4 * eye(n) - diag(ones(n - 1, 1), 1) - diag(ones(n - 1, 1), -1);
L = cell(n, 1);
Ufac = cell(n, 1);
Yrhs = zeros(n, n);
Uint = zeros(n, n);

factor_tic = tic;
Ufac{1} = D;
for i = 2:n
    L{i} = -(Ufac{i - 1} \ eye(n));
    Ufac{i} = D + L{i};
end

Yrhs(:, 1) = F(1, :).';
for i = 2:n
    Yrhs(:, i) = F(i, :).'- L{i} * Yrhs(:, i - 1);
end

Uint(:, n) = Ufac{n} \ Yrhs(:, n);
for i = n-1:-1:1
    Uint(:, i) = Ufac{i} \ (Yrhs(:, i) + Uint(:, i + 1));
end
factor_time = toc(factor_tic);

U = zeros(N + 1, N + 1);
U(2:N, 2:N) = Uint.';

info = struct();
info.h = h;
info.factor_time = factor_time;
end
