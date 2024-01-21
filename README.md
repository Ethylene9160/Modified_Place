# Modified `place` in matlab

![Static matlab](https://img.shields.io/badge/matlab-2022a-green)
![Static cpp](https://img.shields.io/badge/c++-11-blue)
![Static lc](https://img.shields.io/badge/Linear-Control-red)

Switch language to: [中文](sources/Readme_zh.md)

**Project Structure** (major part)

```bash
│  .gitignore
│  LICENSE
│  README.md
├─sources
│	└─Readme_zh.md
├─doc
│	└─doc.md
└─src
	├─cpp
	└─matlab
   		├─myPlace.m
   		├─myOrderedPlace.m
   		├─myRandomPlace.m
   		├─myAdvancedRandomPlace.m
   		└─test_samples.mlx
```

## Introduction

A second project of *SDM364, Multi-variable Control and Applications*, using matlab while basing on cpp to enhance the function `place(A,B,poles)` in matlab.

Before applying this repo, make sure that you're able to access:

> c++ version: standard 11
>
> matlab: R2022a(or replacement)

## Methods

### Controllable SISO Systems

Based on **Bass-Gura Formula**. 

Pass.

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

### Multiple Solutions

For a given piar $(A,B,P)$, letting $T_B$ a non-singular matrix.

Derive:
$$
B' = BT_B
$$
This means, $B'$ is a new column combination of $B$. And let's place pair $(A,B',P)$, assuming what we get is $K'$, and the new state matrix will be: 
$$
A-B'K'\\
=A-(BT_B)\cdot K'\\
=A-B\cdot (T_BK')
$$
Thus, we know, using any non-singular matrix $T_B$, to drive column transform on $B$, and apply row transform of $T_B$ on the obtained $K'$, can we get $K$, That is:
$$
K= T_BK'
$$
That is to say, to find non-singular $T_B$, can get different answers.

However, we're facing 2 problems. The one is that, if we randomly construct a new linear combination of $B$, for $(A,B)$ is controllable, we're easy to see, using the previous methods the first column of $B$ is quite easy to make $(A,b_1)$ controllable, and so we only get a certain series of the solution. So we're going to construct new combination of $B$ in a random but not that random way.

One way is to make simple column exchange. For example, the folloing matrix $T_B$, and $B\cdot T_B$ represents exchanging the first and the second column of matrix $B$:
$$
T_B = \begin{bmatrix}
0~~1~~0\\
1~~0~~0\\
0~~0~~1
\end{bmatrix}
$$
Then we can get a new set of our pole place. Based on this, function `myOrderedPlace(A,B,P,order)` is to make simple combination of matrix $B$, and get new solution set of $K$.

Another way is to make combination in a more complex way. The experimenter use a sercial algorithm to realize such result:
$$
b_i\leftarrow c_1b_1+c_2b_2+...
$$
with $c_i$ either a zero, or a random number in $(0,1)$, The probability of zero is given a certain proportion, and the remaining probability is assigned to a random number between $(0,1)$， to reduce the phenomenon of the apprance of controllable $(A,b_1)$. And finally process it, we get the pole placement feedback matrix $K$

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
>
> 3 For the uncontrollable multi-input matrix, the strategy should be changed to: first, the matrix was split into a controllable part and an uncontrollable part, and then only the controllable part was placed with poles. I'm going to put this part right here in K and I'm going to zero it out, and then I'm going to go back and transform it.(done)

# Update logs

> 24/1/1.  16:33: finished basic SISO with controllable system.
>
> 24/1/1.  19:55: finished SISO with un-controllable system, using controllibility decomposition.
>
> 24/1/2.  15:29: finished basic MIMO sisytems.
>
> 24/1/2.  18:51: fix bugs-modify $B$ to $BT^{-1}$ when passing that matrix to the recursion body.
>
> 24/1/15.21:56: Complete pole configuration under column exchange of input matrix.
>
> 24/1/17.11:10: fix bugs of using SVD to decompose a certain matrix during final exam.
>
> 24/1/18.14:01: New feature - we can place an MIMO system with random column combination of matrix B!
>
> 24/1/18.20:00: Sourcing.



