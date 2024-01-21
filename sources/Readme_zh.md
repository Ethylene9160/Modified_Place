# 覆写`place`函数-matlab

![Static matlab](https://img.shields.io/badge/matlab-2022a-green)
![Static cpp](https://img.shields.io/badge/c++-11-blue)
![Static lc](https://img.shields.io/badge/Linear-Control-red)

将语言切换到：[English](../README.md)

**项目结构**（主要部分）

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

## 简介

这是*SDM364线性系统控制*课程的第二个课程项目，作者将会在matlab中重写place函数，并基于c++为其添加部分功能。

在使用这个仓库前，请保证您满足以下条件：

> c++11及以上标准
>
> matlab 2022a或替代

## 方法

### 可控单输入控制系统设计

基于**Bass-Gura**公式。

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

### 多解性

对于给定的极点配置条件 $(A,B,P)$，设 $T_B$ 为非奇异矩阵。

推导：
$$
B' = BT_B
$$
这里，$B'$ 是 $B$ 的新列组合。使用极点配置条件为 $(A,B',P)$，假设经过极点配置，我们得到的反馈矩阵是 $K'$，那么新的状态矩阵$A$将会是：
$$
A-B'K'\\
=A-(BT_B)\cdot K'\\
=A-B\cdot (T_BK')
$$
因此，我们知道，使用任何非奇异矩阵 $T_B$ 对 $B$ 进行列变换，并对得到的 $K'$ 应用 $T_B$ 的行变换，我们可以得到 $K$，即：
$$
K= T_BK'
$$
也就是说，寻找不同的非奇异的 $T_B$，可以得到不同的答案。

然而，我们面临两个问题。一个是，如果我们随机构造 $B$ 的新线性组合，由于 $(A,B)$ 是可控的，我们容易看出，使用前面的方法容易导致 $B$ 的第一列很容易因为组合了所有的$b_i$，使得 $(A,b_1)$ 已经可控，导致按照实验者所采用的前面的方法进行极点配置的话，仅能得到某一特定系列的解。所以我们要以一种随机但不完全随机的方式构造 $B$ 的新组合。

一种方法是仅仅进行简单的列交换。例如，下面的矩阵 $T_B$，$B\cdot T_B$ 表示交换矩阵 $B$ 的第一列和第二列：
$$
T_B = \begin{bmatrix}
0~~1~~0\\
1~~0~~0\\
0~~0~~1
\end{bmatrix}
$$
然后我们可以得到一组新的极点位置。基于此，函数 `myOrderedPlace(A,B,P,order)` 用于对矩阵 $B$ 进行简单组合，并得到新的 $K$ 解集。

另一种方式是以更复杂的方式进行线性组合。实验者使用特定算法实现这样的结果：
$$
b_i\leftarrow c_1b_1+c_2b_2+...
$$
其中 $c_i$ 是零或 $(0,1)$ 中的随机数，其中零出现的概率具有一定的占比，剩下的概率将分配给$(0,1)$之间的随机数，以减少 $(A,b_1)$ 可控现象的出现。最后对它们进行处理，可以得到极点放置反馈矩阵 $K$。


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

> 1将不可控部分与全零部分进行拼接。（完成）
>
> 2 判断各可控矩阵的维度。例如，如果$A_c$等于$0\times0$，则需要执行其他操作。（完成）
>
> 3 对于不可控多输入矩阵，策略应当更改为，先将矩阵拆分为可控部分和不可控部分，然后仅对可控部分进行极点放置；将这一部分放置好的K后进行补零，再进行回变换。（已完成）

# 更新日志

> 24/1/1.  16:33：完成基本的SISO可控系统。
>
> 24/1/1.  19:55：完成SISO与不可控系统——使用可控性分解。
>
> 24/1/2.  15:29：完成基本的MIMO系统。
>
> 24/1/2.  18:51：修复漏洞-将传递到递归体的原输入矩阵 $B$ 更正为 $BT^{-1}$ 。
>
> 24/1/15.21:56：完成输入矩阵列交换下的极点配置。
>
> 24/1/17.11:10：修复了因SVD分解导致的能控矩阵分解出错的问题。（考 场 科 研）
>
> 24/1/18.14:01：新特性-我们可以用矩阵B的随机线性组合来对MIMO系统进行极点配置!
>
> 24/1/18.20:00：整理。
