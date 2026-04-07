function [x, y, U, info] = poisson_dst_dirichlet(a, b, I, J, f_handle, alpha_handle)

if I < 1 || J < 1 || floor(I) ~= I || floor(J) ~= J
    error('I and J must be positive integers.');
end
if a <= 0 || b <= 0
    error('a and b must be positive.');
end

h = a / (I + 1);
k = b / (J + 1);
x = linspace(0, a, I + 2).';
y = linspace(0, b, J + 2);
xi = x(2:end-1);
yj = y(2:end-1);
[Xi, Yj] = ndgrid(xi, yj);

rhs = f_handle(Xi, Yj);
rhs = double(rhs);
if ~isequal(size(rhs), [I, J])
    error('f_handle must return an array of size I-by-J for interior grid input.');
end

left_bc = reshape(alpha_handle(zeros(1, J), yj), 1, []);
right_bc = reshape(alpha_handle(a * ones(1, J), yj), 1, []);
bottom_bc = reshape(alpha_handle(xi, zeros(I, 1)), [], 1);
top_bc = reshape(alpha_handle(xi, b * ones(I, 1)), [], 1);

rhs(1, :) = rhs(1, :) + left_bc / h^2;
rhs(end, :) = rhs(end, :) + right_bc / h^2;
rhs(:, 1) = rhs(:, 1) + bottom_bc / k^2;
rhs(:, end) = rhs(:, end) + top_bc / k^2;

lambda = 4 / h^2 * sin((1:I)' * pi / (2 * (I + 1))).^2;
mu = 4 / k^2 * sin((1:J) * pi / (2 * (J + 1))).^2;

transform_tic = tic;
V = dst1_fft(dst1_fft(rhs, 1), 2);
W = V ./ (lambda * ones(1, J) + ones(I, 1) * mu);
Uint = (2 / (I + 1)) * (2 / (J + 1)) * dst1_fft(dst1_fft(W, 1), 2);
transform_time = toc(transform_tic);

U = zeros(I + 2, J + 2);
U(1, :) = reshape(alpha_handle(zeros(1, J + 2), y), 1, []);
U(end, :) = reshape(alpha_handle(a * ones(1, J + 2), y), 1, []);
U(:, 1) = reshape(alpha_handle(x, zeros(I + 2, 1)), [], 1);
U(:, end) = reshape(alpha_handle(x, b * ones(I + 2, 1)), [], 1);
U(2:end-1, 2:end-1) = Uint;

info = struct();
info.h = h;
info.k = k;
info.lambda = lambda;
info.mu = mu;
info.rhs = rhs;
info.transform_time = transform_time;
end
