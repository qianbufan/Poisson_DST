function main_project1()

clc;
close all;
output_dir = fullfile(pwd, 'results');
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

%% Task 1
fprintf('[Task 1] Homogeneous Dirichlet on [0,1]^2\n');
f_handle = @(x, y) 2 * pi^2 * sin(pi * x) .* sin(pi * y);
alpha_handle = @(x, y) zeros(size(x));
uexact = @(x, y) sin(pi * x) .* sin(pi * y);
Ns = 10:5:100;
iter_opts = struct('tol', 1e-8, 'maxIter', 50000);
rows = {};

for idx = 1:numel(Ns)
    N = Ns(idx);
    I = N - 1;
    J = N - 1;

    tic;
    [x, y, U_dst] = poisson_dst_dirichlet(1, 1, I, J, f_handle, alpha_handle);
    t_dst = toc;
    [err_inf, err_rel] = compute_error(x, y, U_dst, uexact);
    rows(end + 1, :) = {'DST', N, t_dst, NaN, err_inf, err_rel}; 

    tic;
    [x, y, U_blk] = poisson_block_lu_square(N, f_handle);
    t_blk = toc;
    [err_inf, err_rel] = compute_error(x, y, U_blk, uexact);
    rows(end + 1, :) = {'BlockLU', N, t_blk, NaN, err_inf, err_rel}; 

    tic;
    [x, y, U_jac, info_jac] = poisson_jacobi_dirichlet(1, 1, I, J, f_handle, alpha_handle, iter_opts);
    t_jac = toc;
    [err_inf, err_rel] = compute_error(x, y, U_jac, uexact);
    rows(end + 1, :) = {'Jacobi', N, t_jac, info_jac.iter, err_inf, err_rel}; 

    tic;
    [x, y, U_gs, info_gs] = poisson_gs_dirichlet(1, 1, I, J, f_handle, alpha_handle, iter_opts);
    t_gs = toc;
    [err_inf, err_rel] = compute_error(x, y, U_gs, uexact);
    rows(end + 1, :) = {'Gauss-Seidel', N, t_gs, info_gs.iter, err_inf, err_rel}; 
end

T1 = cell2table(rows, 'VariableNames', {'Method', 'N', 'TimeSeconds', 'Iterations', 'MaxError', 'RelativeFroError'});
disp(T1);
writetable(T1, fullfile(output_dir, 'task1_comparison.csv'));
plot_task1(Ns, T1, output_dir);

%% Task 2
fprintf('\n[Task 2] Non-homogeneous Dirichlet on [0,1]^2\n');
a = 1;
b = 1;
I = 127;
J = 127;
uexact2 = @(x, y) exp(x + y);
f_handle2 = @(x, y) -2 * exp(x + y);
alpha_handle2 = @(x, y) exp(x + y);
[x2, y2, U2] = poisson_dst_dirichlet(a, b, I, J, f_handle2, alpha_handle2);
[err_inf2, err_rel2, Ue2] = compute_error(x2, y2, U2, uexact2);
fprintf('Task 2 max error        : %.6e\n', err_inf2);
fprintf('Task 2 relative fro err : %.6e\n', err_rel2);
plot_surface_triplet(x2, y2, U2, Ue2, 'task2_nonhomogeneous', output_dir);

%% Task 3
fprintf('\n[Task 3] General rectangle [0,a] x [0,b]\n');
a = 2.0;
b = 1.5;
I = 100;
J =75;
uexact3 = @(x, y) sin(pi * x / a) .* sin(2 * pi * y / b);
f_handle3 = @(x, y) ((pi / a)^2 + (2 * pi / b)^2) .* sin(pi * x / a) .* sin(2 * pi * y / b);
alpha_handle3 = @(x, y) zeros(size(x));
[x3, y3, U3] = poisson_dst_dirichlet(a, b, I, J, f_handle3, alpha_handle3);
[err_inf3, err_rel3, Ue3] = compute_error(x3, y3, U3, uexact3);
fprintf('Task 3 max error        : %.6e\n', err_inf3);
fprintf('Task 3 relative fro err : %.6e\n', err_rel3);
plot_surface_triplet(x3, y3, U3, Ue3, 'task3_rectangle', output_dir);

fprintf('\nAll tables and figures were saved to: %s\n', output_dir);
end

function [err_inf, err_rel, Ue] = compute_error(x, y, U, uexact)
[X, Y] = ndgrid(x, y);
Ue = double(uexact(X, Y));
D = U - Ue;
err_inf = max(abs(D(:)));
inner_num = D(2:end-1, 2:end-1);
inner_exact = Ue(2:end-1, 2:end-1);
err_rel = norm(inner_num, 'fro') / max(norm(inner_exact, 'fro'), eps);
end

function plot_task1(Ns, T, output_dir)
method_names = unique(T.Method, 'stable');
figure('Color', 'w', 'Name', 'Task 1 comparison');

subplot(1, 2, 1);
hold on;
for i = 1:numel(method_names)
    mask = strcmp(T.Method, method_names{i});
    semilogy(T.N(mask), T.TimeSeconds(mask), '-', 'LineWidth', 1.5, 'DisplayName', method_names{i});
end
xlabel('N');
ylabel('Time (s)');
title('Runtime comparison');
grid on;
legend('Location', 'northwest');

subplot(1, 2, 2);
hold on;

for i = 1:numel(method_names)
    mask = strcmp(T.Method, method_names{i});
    semilogy(T.N(mask), T.MaxError(mask), '-', 'LineWidth', 1.5, 'DisplayName', method_names{i});
end
xlabel('N');
ylabel('Max error');
title('Accuracy comparison');
grid on;
legend('Location', 'northeast');

ax_main = gca;
ax_inset = axes('Position', [0.68, 0.35, 0.22, 0.35]);  

mask_dst = strcmp(T.Method, 'DST');
semilogy(ax_inset, T.N(mask_dst), T.MaxError(mask_dst), '-', 'LineWidth', 1.5, 'Color', [0, 0.4470, 0.7410]);
title(ax_inset, 'DST detail');
grid(ax_inset, 'on');
box(ax_inset, 'on');

set(ax_inset, 'FontSize', 8);
set(ax_inset, 'Color', [1 1 0.95]);  

saveas(gcf, fullfile(output_dir, 'task1_comparison.png'));
end

function plot_surface_triplet(x, y, U_num, U_exact, stem_name, output_dir)
[X, Y] = ndgrid(x, y);
Err = abs(U_num - U_exact);
figure('Color', 'w', 'Name', stem_name, 'Position', [100, 100, 1200, 400]);
subplot(1, 3, 1);
surf(X, Y, U_num, 'EdgeColor', 'none');
view(45, 30);
xlabel('x'); ylabel('y'); zlabel('u_h');
title('Numerical solution'); colorbar;

subplot(1, 3, 2);
surf(X, Y, U_exact, 'EdgeColor', 'none');
view(45, 30);
xlabel('x'); ylabel('y'); zlabel('u');
title('Exact solution'); colorbar;

subplot(1, 3, 3);
surf(X, Y, Err, 'EdgeColor', 'none');
view(45, 30);
xlabel('x'); ylabel('y'); zlabel('|e|');
title('Absolute error'); colorbar;
saveas(gcf, fullfile(output_dir, [stem_name, '.png']));
end
