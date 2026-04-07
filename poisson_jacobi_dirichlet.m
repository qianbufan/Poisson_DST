function [x, y, U, info] = poisson_jacobi_dirichlet(a, b, I, J, f_handle, alpha_handle, opts)

if nargin < 7
    opts = struct();
end
if ~isfield(opts, 'tol'), opts.tol = 1e-8; end
if ~isfield(opts, 'maxIter'), opts.maxIter = 50000; end
if ~isfield(opts, 'U0'), opts.U0 = []; end

h = a / (I + 1);
k = b / (J + 1);
den = 2 / h^2 + 2 / k^2;
x = linspace(0, a, I + 2).';
y = linspace(0, b, J + 2);
xi = x(2:end-1);
yj = y(2:end-1);
[Xi, Yj] = ndgrid(xi, yj);
F = double(f_handle(Xi, Yj));

U = zeros(I + 2, J + 2);
U(1, :) = reshape(alpha_handle(zeros(1, J + 2), y), 1, []);
U(end, :) = reshape(alpha_handle(a * ones(1, J + 2), y), 1, []);
U(:, 1) = reshape(alpha_handle(x, zeros(I + 2, 1)), [], 1);
U(:, end) = reshape(alpha_handle(x, b * ones(I + 2, 1)), [], 1);
if ~isempty(opts.U0)
    U(2:end-1, 2:end-1) = opts.U0;
end

V = U;
iter_tic = tic;
for iter = 1:opts.maxIter
    for i = 2:I+1
        for j = 2:J+1
            V(i, j) = ((U(i - 1, j) + U(i + 1, j)) / h^2 + (U(i, j - 1) + U(i, j + 1)) / k^2 + F(i - 1, j - 1)) / den;
        end
    end
    tmp = V(2:end-1, 2:end-1) - U(2:end-1, 2:end-1);
    diff_now = max(abs(tmp(:)));
    U = V;
    if diff_now < opts.tol
        break;
    end
end
elapsed = toc(iter_tic);

info = struct();
info.iter = iter;
info.tol = opts.tol;
info.maxIter = opts.maxIter;
info.diff = diff_now;
info.elapsed = elapsed;
info.converged = diff_now < opts.tol;
end
