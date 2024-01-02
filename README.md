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

## Design MIMO Model

For higher dimension matrix:
$$
B=\begin{bmatrix}
b_1~b_2~b_3...b_{c_b}
\end{bmatrix}\\
$$
Let's start from the first column of $B$. If $\begin{bmatrix}
A~|~b_1
\end{bmatrix}$ is controllable, then we only need to calculate $b_1$, and we can get the first row of $K$.

If $\begin{bmatrix}
A~|~b_1
\end{bmatrix}$ is not controllable, in the [last chapter](###Uncontrollable SISO Systems), we have discussed the decomposition of the controlling matrixes. So if we apply decomposition on $\begin{bmatrix}
A~|~b_1
\end{bmatrix}$, then we get the new controlling matrixes corresponding to $T$:
$$
\begin{bmatrix}
A~|~b_1
\end{bmatrix}
\rightarrow
\begin{bmatrix}
A_{11}~A_{12}~|~b_{11}\\
~0~~~A_{22}~|~~0~
\end{bmatrix}
$$
The size of $T$ is $n\times n$, so we can apply it on matrix $B$, whose dimension is $n\times c_b$. The new matrix $B$ can be written as:
$$
\hat B = T^{-1}B=\begin{bmatrix}
b_{11}~b_{12}~b_{13}...\\
~0~~~b_{22}~b_{23}...
\end{bmatrix}
$$
So finally, the comtrolling pair can be written as:
$$
\begin{bmatrix}
A~|~B
\end{bmatrix}
\rightarrow
\begin{bmatrix}
A_{11}~A_{12}~|~b_{11}~B_1\\
~0~~~A_{22}~|~~0~~B_2
\end{bmatrix}
$$
Assume the rank of the controllable part $A_{11}$ is $r$, then the vector $K$(corresponding to this row) can be get using $A_{11},~b_{11}$ and the first $r$ columns of target poles.

Then, it's a good way to CONTINUE using pair $\begin{bmatrix}A_{22}~|~B_2\end{bmatrix}$, to go on with this. Assume we now get a new $K_{new}$. We can simply add it to the next block position of $K$:
$$
K\leftarrow \begin{bmatrix}
K~~~~~~0~~~\\
~0~~~~K_{new}
\end{bmatrix}
$$
And if $\begin{bmatrix}A~|~B_2\end{bmatrix}$ is **controllable** and we successfully get the $K$, we can simply append zeros on the buttom of $K$ until the row of $K$ is the same value of $c_b$(colums of matrix $B$) and finish calculating $K$.
$$
K\leftarrow \begin{bmatrix}
K~~~~~~0~~~\\
~0~~~~K_{new}\\
~0~\\
...
\end{bmatrix}
$$
And if  $\begin{bmatrix}A_{22}~|~B_2\end{bmatrix}$ is **uncontrollable**, we can also do decomposition and use piar  $\begin{bmatrix}A_{33}~|~B_3\end{bmatrix}$, till the we get the final $K$.

The author use recursion to realize this.

## preponderance

### Pole placement for uncontrollable performance control matrices

The built-in  function `place` in matlab does not support the input of uncontrollable controllability matrix pairs. The modified function realizes the pole placement of uncontrollable controllability matrix by decomposing the controllability matrix.

### No more restrictions on geometric multiplicity

The built-in  function `place` in matlab limits the geometric multiplicity of the poles. If the geometric multiplicity of the pole matrix is greater than the geometric multiplicity of the controllability matrix, the function will not work. By using the traditional formula decomposition method, the overwritten function does not have the condition restriction of geometric multiplicity.

## shortcomes

### Lower speed

Most operations are operated in matlab rather than lower-level computitional languages like `c` or `cpp`, which restricts the speed of calculating process.

## TODO

> 1 Joint uncontrollable part together with all-zero part.(fixed)
>
> 2 Judge dimension of each controllable part. For example, if $A_c$ is $0\times 0$, we need to do other operations.(fixed)

# Update logs

> 24/1/1.16:33: finished basic SISO with controllable system.
>
> 24/1/1.19:55: finished SISO with un-controllable system, using controllibility decomposition.
>
> 24/1/2.15:29: finished basic MIMO sisytems.
>
> 24/1/2.18:51: fix bugs-modify $B$ to $BT^{-1}$ when passing that matrix to the recursion body.

$$
\begin{bmatrix}
A~|~b_1~b_2
\end{bmatrix}

\begin{bmatrix}
A_{11}~A_{12}~|~b_{11}~b_{12}\\
~0~~~A_{22}~|~0~~~b_{22}
\end{bmatrix}
$$

```matlab
function K = myPoleMIMOPlacementRecursion(A,B,P)
    if size(B,2) == 1
        disp("final task!")
        K = myPoleSISOPlacement(A,B,P)
    else
        disp("current input B is:")
        disp(B)
        if size(B,2) == 1
            Bi = B;
        else
            Bi = B(:,1:1);
        end
        [T, A_c, B_c, A_uc, B_uc] = controllabilityDecomposition(A,Bi);
        disp("size of A_c is:")
        r = size(A_c, 1)
        if r ~= 0
            K_up = myPoleSISOPlacement(A_c, B_c, P(1:r))
        else
            K_up=[]
        end
        A = A_uc;
        %B = B_uc(:,2:end)
        B = B(r+1:end,2:end);
        P = P(r+1:end);
        r;
        s = size(A,1);
        if size(A, 1) == 0
            str = "over!"
            %k_up_size= size(K_up,2)
            %K = [K_up;zeros(1,size(K_up,2))]
            K_up = [K_up;zeros(1,size(K_up,2))];
            K = K_up*inv(T)
            
        else
            
            str = "continue!"
            K_down = myPoleMIMOPlacementRecursion(A,B,P)
            if r == 0
                K = [zeros(1, size(K_down,2));
                    K_down];
            else
                K = [K_up, zeros(size(K_up,1), size(K_down,2));
                    zeros(size(K_down,1), size(K_up, 2)), K_down]
            end
            K=K*inv(T)
        end
    end
end
```

