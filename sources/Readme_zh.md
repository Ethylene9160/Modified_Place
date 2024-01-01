# 覆写`place`函数-matlab

![Static matlab](https://img.shields.io/badge/matlab-2022a-green)

![Static cpp](https://img.shields.io/badge/c++-11-blue)

![Static lc](https://img.shields.io/badge/Linear-Control-red)

将语言切换到：[English](../README.md)

**项目结构**

```bash
│  .gitignore
│  LICENSE
│  README.md
├─sources
│  ├─Readme_zh.md
└─src
   ├─cpp
   └─matlab
```

## 简介

这是SDM364线性系统控制课程的第二个课程项目，作者将会在matlab中重写place函数，并基于c++为其添加部分功能。

在使用这个仓库前，请保证您满足以下条件：

> c++11及以上标准
>
> matlab 2022a或替代

## 方法

### 可控SISO系统

略。

### 不可控SISO系统

将系统分解为可控部分和不可控部分。例如：
$$
A=\begin{bmatrix}
A_{11}~A_{12}\\~0~~~A_{22}
\end{bmatrix}\\
B=\begin{bmatrix}
B_{11}\\0
\end{bmatrix}
$$
这时候，需要寻找转换矩阵T。

在对能控矩阵进行分解的时候，我们知道，$A=T^{-1}\Lambda T$，其中$\Lambda$

是一个对角矩阵。若矩阵不可控，显然不存在这样的$\Lambda$，但是我们可以通过某些手段将其分解为约当型下三角矩阵。

计算出$M_c$，显然它不是满秩的（假设秩为$r<n$），它不可逆。因此，我们可以提取$M_c$中线性无关的$r$列，然后将其填充到$T$的前$r$列。

显然，我们想要求出转换矩阵的逆变换，必然要求转换矩阵满秩。因此我们需要在$T$的后$n-r$列填充上让转换矩阵$T$满秩的列。

这样的列怎么找呢？显然，左零空间反映了矩阵列的线性相关情况，因而我们求出$T$的左零空间，将左零空间补充到这$n-r$列，就能实现分解了。

```matlab
function [T, A_c, B_c, A_uc, B_uc] = controllabilityDecomposition(A, B)
    % Calculate the controllability matrix
    C = ctrb(A, B);
    
    % Determine the rank of the controllability matrix
    rankC = rank(C);
    
    % Select rankC independent columns from C
    [U, S, V] = svd(C, 'econ');
    independent_cols = V(:, 1:rankC);
    
    % Ensure we have a full set of basis vectors for the state space
    % If system is not fully controllable, complete the basis.
    if rankC < size(A,1)
        % Add additional vectors to span the entire space
        % Find the null space of the controllability matrix
        null_vectors = null(C');
        
        % Construct T by combining the independent columns and the null space
        T = [independent_cols, null_vectors];
    else
        % System is fully controllable
        T = independent_cols;
    end
    
    % Verify that T is full rank (non-singular)
    if rank(T) < size(A,1)
        error('Transformation matrix T is singular. The system may not be controllable.');
    end
    
    % Transform the system matrices
    T_inv = inv(T);
    A_tilde = T_inv * A * T;
    B_tilde = T_inv * B;
    
    % Extract the submatrices for the controllable part
    A_c = A_tilde(1:rankC, 1:rankC);
    B_c = B_tilde(1:rankC, :);
    
    % Extract the submatrices for the uncontrollable part, if any
    A_uc = A_tilde(rankC+1:end, rankC+1:end);
    B_uc = B_tilde(rankC+1:end, :);
end
```

