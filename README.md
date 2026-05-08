# Poisson方程五点差分格式的快速算法与数值比较实验

本项目按照 `Project-1.pdf` 的要求实现了二维 Poisson 方程五点差分格式的快速求解与数值实验，并参考 `FFTPoissonCode.pdf` 中给出的 DST/BlockLU/Jacobi/Gauss-Seidel 思路。

## 1. 已实现内容

1. **在单位正方形 `[0,1] x [0,1]` 上，齐次 Dirichlet 边界条件下的快速 DST 求解器**，并与 Block LU、Jacobi、Gauss-Seidel 方法做时间和误差比较。
2. **扩展到非齐次 Dirichlet 边界条件**。
3. **扩展到一般矩形区域 `[0,a] x [0,b]`**。

## 2. 文件说明

```text
project1_matlab/
├── main_project1.m              % 主程序：一键运行全部实验
├── poisson_dst_dirichlet.m      % 通用 DST 快速求解器，支持矩形区域 + 非齐次 Dirichlet
├── poisson_block_lu_square.m    % [0,1]^2、零边界条件下的 Block LU 对比程序
├── poisson_jacobi_dirichlet.m   % Jacobi 迭代法
├── poisson_gs_dirichlet.m       % Gauss-Seidel 迭代法
├── dst1_fft.m                   % 基于 FFT 的 DST-I 实现
└── README.md                    % 使用说明
```

## 3. 运行方法

将整个文件夹加入 MATLAB 工作路径后，在命令行运行：

```matlab
main_project1
```

程序会自动完成三部分实验，并在当前目录下生成：

```text
results/
├── task1_comparison.csv
├── task1_comparison.png
├── task1_dst_error.png
├── task2_nonhomogeneous.png
└── task3_rectangle.png
```

## 4. 数值模型说明

求解问题为：

```math
-\Delta u = f \quad \text{in } \Omega,
\qquad
u = \alpha \quad \text{on } \partial\Omega.
```

离散区域取为：

```math
\Omega = [0,a] \times [0,b],
\qquad
h = \frac{a}{I+1},
\qquad
k = \frac{b}{J+1}.
```

内部未知量满足五点差分格式：

```math
\left(\frac{2}{h^2} + \frac{2}{k^2}\right)u_{ij}
- \frac{u_{i-1,j}+u_{i+1,j}}{h^2}
- \frac{u_{i,j-1}+u_{i,j+1}}{k^2}
= f_{ij}.
```

对于 Dirichlet 边界，边界项被直接并入右端项，因此 `poisson_dst_dirichlet.m` 可统一处理：

- 齐次边界；
- 非齐次边界；
- 一般矩形区域。

## 5. 主程序中的测试问题

### Task 1: `[0,1]^2` 上齐次 Dirichlet

精确解：

```matlab
u(x,y) = sin(pi*x).*sin(pi*y)
```

对应右端项：

```matlab
f(x,y) = 2*pi^2*sin(pi*x).*sin(pi*y)
```

边界条件：

```matlab
alpha(x,y) = 0
```

该问题用于：
- 验证快速 DST 求解器的正确性；
- 比较 DST、BlockLU、Jacobi、Gauss-Seidel 的运行时间与误差。
<img width="875" height="656" alt="task1_comparison" src="https://github.com/user-attachments/assets/c96a7498-df7f-47c6-b3b2-bec323bf4740" />

### Task 2: `[0,1]^2` 上非齐次 Dirichlet

精确解：

```matlab
u(x,y) = exp(x+y)
```

对应右端项：

```matlab
f(x,y) = -2*exp(x+y)
```

边界条件：

```matlab
alpha(x,y) = exp(x+y)
```

该问题用于验证非齐次边界处理是否正确。
<img width="1875" height="625" alt="task2_nonhomogeneous" src="https://github.com/user-attachments/assets/b2ade596-a5da-46fd-8bf4-f83d8b717ba5" />


### Task 3: 一般矩形 `[0,a] x [0,b]`

主程序中取：

```matlab
a = 2.0;
b = 1.5;
```

精确解：

```matlab
u(x,y) = sin(pi*x/a).*sin(2*pi*y/b)
```

对应右端项：

```matlab
f(x,y) = ((pi/a)^2 + (2*pi/b)^2) .* sin(pi*x/a).*sin(2*pi*y/b)
```

边界条件仍为零边界，用来验证程序已不依赖单位正方形假设。
<img width="1875" height="625" alt="task3_rectangle" src="https://github.com/user-attachments/assets/4fe0f0e8-6ff1-4890-82c6-27ba692023d1" />


## 6. 关于 DST 实现

`dst1_fft.m` 实现的是 **DST-I**，采用 FFT 构造奇延拓完成。

若一维向量长度为 `n`，则其逆变换关系为：

```matlab
x = 2/(n+1) * dst1_fft(y)
```

二维求解中，对两个方向分别施加 DST-I，并利用离散 Laplace 算子的特征值：

```math
\lambda_i = \frac{4}{h^2} \sin^2\left(\frac{i\pi}{2(I+1)}\right),
\qquad
\mu_j = \frac{4}{k^2} \sin^2\left(\frac{j\pi}{2(J+1)}\right).
```

## 7. 输出指标

主程序输出并保存：

- `TimeSeconds`：求解耗时；
- `Iterations`：迭代法步数，直接法记为 `NaN`；
- `MaxError`：全网格最大绝对误差；
- `RelativeFroError`：内部点相对 Frobenius 误差。

