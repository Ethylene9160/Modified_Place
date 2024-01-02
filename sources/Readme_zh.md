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
│	└─Readme_zh.md
└─src
	├─cpp
	└─matlab
   		├─myPoleSISOPlacement.m
   		├─test_samples.mlx
   		└─untitled.mlx
```

## 简介

这是*SDM364线性系统控制*课程的第二个课程项目，作者将会在matlab中重写place函数，并基于c++为其添加部分功能。

在使用这个仓库前，请保证您满足以下条件：

> c++11及以上标准
>
> matlab 2022a或替代

## 方法

### 可控单输入控制系统设计

略。

### 不可控单输入控制系统设计

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

是一个对角矩阵。若矩阵不可控，显然不存在这样的纯对角线矩阵$\Lambda$，但是我们可以通过某些手段将其分解为约当型下三角矩阵。

计算出$M_c$，显然它不是满秩的（假设秩为$r<n$），它不可逆。因此，我们可以提取$M_c$中线性无关的$r$列，然后将其填充到$T$的前$r$列。

显然，我们想要求出转换矩阵的逆变换，必然要求转换矩阵满秩。因此我们需要在$T$的后$n-r$列填充上让转换矩阵$T$满秩的列。补充了零空间的转换矩阵$T$可以进行约当分解。

这样的列怎么找呢？显然，左零空间反映了矩阵列的线性相关情况，因而我们求出$T$的左零空间，将左零空间补充到这$n-r$列，就能实现分解了。

```matlab
function [T, A_c, B_c, A_uc, B_uc] = controllabilityDecomposition(A, B)
    % 能控矩阵
    C = ctrb(A, B);
    
    % 确定可控性矩阵的秩
    rankC = rank(C);
    
    % 从C中选择rankC独立列
    [U, S, V] = svd(C, 'econ');
    independent_cols = V(:, 1:rankC);
    
    % 确保我们有状态空间的一整套基向量。
    % 如果系统不完全可控，则需要手动计算转换矩阵。
    if rankC < size(A,1)
        % 添加额外的向量来张成整个空间
        % 找到可控性矩阵的左零空间
        null_vectors = null(C');
        
        % 通过合并独立列和左零空间来构造T
        T = [independent_cols, null_vectors];
    else
        % 若系统本来可控
        T = independent_cols;
    end
    
    % 检测T是否奇异（大概率不会）
    if rank(T) < size(A,1)
        error('Transformation matrix T is singular. The system may not be controllable.');
    end
    
    % 变换系统矩阵
    T_inv = inv(T);
    A_tilde = T_inv * A * T;
    B_tilde = T_inv * B;
    
    % 提取可控部分的子矩阵
    A_c = A_tilde(1:rankC, 1:rankC);
    B_c = B_tilde(1:rankC, :);
    
    % 如果有，提取不可控部分的子矩阵
    A_uc = A_tilde(rankC+1:end, rankC+1:end);
    B_uc = B_tilde(rankC+1:end, :);
end
```

完成分解后，我们使用$A_{11}~B_{11}$作为能控的部分，这样能够求出一个维度是$r$的$K_{ccf}$。然后我们在$K_{ccf}$后补零直到$K_{ccf}$的维度变为$n$。将它当作普通的$K_{ccf}$与转换矩阵的逆进行右乘，就能得到：
$$
K = K_{CCF}T^{-1}
$$

## 多输入系统控制器设计

对于高维矩阵
$$
B=\begin{bmatrix}
b_1~b_2~b_3...
\end{bmatrix}\\
$$
如果：$\begin{bmatrix}
A~|~b_1
\end{bmatrix}$能控，那么我们只需要计算出$b_1$，我们可以得到$K$矩阵的第一行。

如果$\begin{bmatrix}A~ | ~ b_1\end{bmatrix}$是不可控的，[上一节](###不可控单输入控制系统设计)我们提到过，我们可以对它采取能控矩阵分解的方式，将其拆分为可控部分和不可控部分。因而对 $\begin{bmatrix}
A~|~b_1
\end{bmatrix}$进行能控分解, 以获取转换矩阵 $T$，并能够得出:
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
$T$的大小是$n\times n$，因此我们可以将其应左乘在矩阵$B$（维数是$n\time c_b$）上。新矩阵$B$可以写成:
$$
\hat B = T^{-1}B=\begin{bmatrix}
b_{11}~b_{12}~b_{13}...\\
~0~~~b_{22}~b_{23}...
\end{bmatrix}
$$
所以最后，能控矩阵对可以写成:
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
假设可控部分$A_{11}$的秩为$r$，则可以用$A_{11}、~b_{11}$和目标极点的前$r$列得到对应这一行的$K$。

然后，继续使用能控矩阵对 $\begin{bmatrix} A_{22} ~|~B_2\end{bmatrix}$。假设我们现在得到了一个新的$K_{new}$。我们可以简单地将其添加到$K$的下一个约当块位置:
$$
K\leftarrow \begin{bmatrix}
K~~~~~~0~~~\\
~0~~~~K_{new}
\end{bmatrix}
$$
如果$\begin{bmatrix}A~|~B_2\end{bmatrix}$是**可控的**，这使得我们成功地得到了$K$，我们可以简单地在$K$的底部添加0，直到$K$的行数$c_b$(矩阵$B$的列)的值相同，结束求取K的过程。
$$
K\leftarrow \begin{bmatrix}
K~~~~~~0~~~\\
~0~~~~K_{new}\\
~0~\\
...
\end{bmatrix}
$$
如果$\begin{bmatrix}A~|~B_2\end{bmatrix}$是**不可控的**，我们则继续时候用能控矩阵对 $\begin{bmatrix}A_{33}~|~B_3\end{bmatrix}$进行分解，直到我们得到最终的$K$。

作者使用递归来实现这一点。


## 优势

经过重写的函数，相对于原有的`place`函数，具有以下优势：

### 非可控性能控矩阵的极点放置

matlab中内置的`place`函数不支持不可控能控矩阵对的输入，经过覆写的函数通过对能控矩阵进行分解，实现了非可控的能控矩阵极点放置。

### 不再对几何重数进行限制

matlab中内置的`place`函数限制了极点的几何重数。如果极点矩阵的几何重数大于能控矩阵的几何重数，该函数将不会工作。经过覆写的函数通过使用传统公式分解法，没有几何重数的条件限制。

## 劣势

### 速度性能的下降

大多数运算都是在matlab中进行的，而不是像`c`或`cpp`这样的底层计算语言，这限制了计算速度。

## 待办

> 1将不可控部分与全零部分进行拼接。
>
> 2 判断各可控矩阵的维度。例如，如果$A_c$等于$0\times0$，则需要执行其他操作。

# 更新日志

> 24/1/1.16:33：完成基本的SISO可控系统。
>
> 24/1/1.19:55：完成SISO与不可控系统——使用可控性分解。
>
> 24/1/2.15:29：完成基本的MIMO系统。
