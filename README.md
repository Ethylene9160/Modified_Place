# Modified `place` in matlab

![Static matlab](https://img.shields.io/badge/matlab-2022a-green)
![Static cpp](https://img.shields.io/badge/c++-11-blue)
![Static lc](https://img.shields.io/badge/Linear-Control-red)

Switch language to: [中文](sources/Readme_zh.md)

**Project Structure**

```bash
│  .gitignore
│  LICENSE
│  README.md
├─sources
│	└─Readme_zh.md
└─src
	├─cpp
	└─matlab
   		├─myPoleSISOPlacement.m
   		├─test_samples.mlx
   		└─untitled.mlx
```

## Introduction

A second project of *SDM364, Multi-variable Control and Applications*, using matlab while basing on cpp to enhance the function `place(A,B,poles)` in matlab.

Before applying this repo, make sure that you're able to access:

> c++ version: standard 11
>
> matlab: R2022a(or replacement)

## Methods

### Controllable SISO Systems

pass.

### Uncontrollable SISO Systems

The system is decomposed into controllable and uncontrollable parts. For example:
$$
A=\begin{bmatrix}
A_{11}~A_{12}\\~0~~~A_{22}
\end{bmatrix}\\
B=\begin{bmatrix}
B_{11}\\0
\end{bmatrix}
$$
In this case, we need to find the transformation matrix $T$.

When we decompose the controllability matrix, we know that, $A=T^{-1}\Lambda T$, with $\Lambda$ a diagnal matrix. If the matrix is not controllable, obviously there is no such matrix $\Lambda$, But we can decompose it into a Jordan lower triangular matrix by some means.

Calculate $M_c$, which is not full rank(assume $r<n$) that is uninvertable. Thus, we can extract the linearly independent $r$ columns in $M_c$,  then padded to the first $r$ column of $T$.

How to find a column like this? Obviously, the left null space reflects the linear dependence of the columns of the matrix, so we can find the left null space of $T$, add the left null space to these $n-r$ columns, and we can realize the decomposition.

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
    
    % Verify that T is full rank (mostly it is full rank)
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

## preponderance

### Pole placement for uncontrollable performance control matrices

The built-in  function `place` in matlab does not support the input of uncontrollable controllability matrix pairs. The modified function realizes the pole placement of uncontrollable controllability matrix by decomposing the controllability matrix.

### No more restrictions on geometric multiplicity

The built-in  function `place` in matlab limits the geometric multiplicity of the poles. If the geometric multiplicity of the pole matrix is greater than the geometric multiplicity of the controllability matrix, the function will not work. By using the traditional formula decomposition method, the overwritten function does not have the condition restriction of geometric multiplicity.

## shortcomes

# Lower speed

Most operations are operated in matlab rather than lower-level computitional languages like `c` or `cpp`, which restricts the speed of calculating process.

# Update logs

> 24/1/1.16:33: finished basic SISO with controllable system.
>
> 24/1/1.19:55: finished SISO with un-controllable system, using controllibility decomposition.

$$
\begin{bmatrix}
A~|~b_1~b_2
\end{bmatrix}

\begin{bmatrix}
A_{11}~A_{12}~|~b_{11}~b_{12}\\
~0~~~A_{22}~|~0~~~b_{22}
\end{bmatrix}
$$

