---
title: "线性模型的理论与求解"
date: 2018-04-20
toc: true
categories:
  - 笔记
tags:
  - 多元统计
  - 回归分析
---

本文总结了线性模型的主要知识点，分别为参数估计，包括最小二乘估计和极大似然估计，区间估计，假设检验。此外，针对每个内容，本文还给出了相应的`R`软件求解算法，并做了相应的模拟。

设在线性模型中：$y\in \mathbb{R}^n$，$X\in \mathbb{R}^{n\times p}$，$n$ 表示观测数，$p$ 表示变量数。$\beta \in \mathbb{R}^n$ 表示回归系数。$\epsilon \sim N\left( 0,\sigma ^2I_n \right) $，是为独立同分布的高斯随机变量。


## 1 估计量的表达式
### 1.1 最小二乘估计量
- **有截距项的参数估计：**

首先令带有截距项的线性模型表达式为：

$$
y=1\beta _0+X\beta+\epsilon 
$$

根据最小二乘思想可知目标函数为 $f\left( \beta _0,\beta  \right) =\lVert y-1\beta _0-X\beta  \rVert _{2}^{2}
$，所以将 $f\left( \beta _0,\beta  \right) $ 展开有：

$$
\begin{split}
f\left( \beta _0,\beta \right) =&y^{\mathrm{T}}y+\beta _01^{\mathrm{T}}1\beta _0+2\beta ^{\mathrm{T}}X^{\mathrm{T}}1\beta _0
\newline
&+\beta ^{\mathrm{T}}X^{\mathrm{T}}X\beta -2\beta _01^{\mathrm{T}}y-2y^{\mathrm{T}}X\beta
\end{split}
$$

然后分别对 $\beta_0$，$\beta$ 求偏导并令为零：

$$
\begin{split}
\frac{\partial f}{\beta _0}&=2\cdot 1^{\mathrm{T}}1\beta _0+2\beta ^{\mathrm{T}}X ^{\mathrm{T}}1-2y ^{\mathrm{T}}1=0
\newline
\frac{\partial f}{\beta }&=2X^{\mathrm{T}}1\beta _0+2X^{\mathrm{T}}X\beta -2X ^{\mathrm{T}}y=0
\end{split}
$$

由上式第1式求得：$\hat{\beta}_0=\left( 1^{\mathrm{T}}1 \right) ^{-1}1^{\mathrm{T}}y-\left( 1^{\mathrm{T}}1 \right) ^{-1}1^{\mathrm{T}}X\beta =\bar{y}-\bar{x}^{\mathrm{T}}\hat{\beta} $，其中 $\bar{X}=\left( \bar{x}_1,\cdots ,\bar{x}_p \right) ^{\mathrm{T}}$ 表示由 $X$ 的每一列的均值组成的列向量。然后将求得的结果带入到第2式中，有：

$$
X^{\mathrm{T}}1\left( 1^{\mathrm{T}}1 \right) ^{-1}1^{\mathrm{T}}y-X^{\mathrm{T}}1\left( 1^{\mathrm{T}}1 \right) ^{-1}1^{\mathrm{T}}X\beta +X^{\mathrm{T}}X\beta -X^{\mathrm{T}}y=0
$$

整理可得：

$$
-X^{\mathrm{T}}\left[ I-1\left( 1^{\mathrm{T}}1 \right) ^{-1}1^{\mathrm{T}} \right] y+X^{\mathrm{T}}\left[ I-1\left( 1^{\mathrm{T}}1 \right) ^{-1}1^{\mathrm{T}} \right] X\beta =0
$$

注意到矩阵 $I-1\left( 1^{\mathrm{T}}1 \right) ^{-1}1^{\mathrm{T}}$ 是对称幂等的，并且 $\left( I-1\left( 1^{\mathrm{T}}1 \right) ^{-1}1^{\mathrm{T}} \right) X$ 可以看成是 $X$ 的每个元素减去所在列的均值得到的新矩阵（就是每一列元素进行中心化）。则记

$$
X_c=\left[ I-1\left( 1^{\mathrm{T}}1 \right) ^{-1}1^{\mathrm{T}} \right] X
$$

所以可以得到

$$
X_{c}^{\mathrm{T}}X_c\beta =X_{c}^{\mathrm{T}}y
$$

故而得到有截距的最小二乘参数估计为：

$$
\begin{split}
\hat{\beta}_0&=\bar{y}-\bar{X}^{\mathrm{T}}\hat{\beta}
\newline
\hat{\beta}&=\left( X _{c}^{\mathrm{T}}X_c \right) ^{-1}X _{c}^{\mathrm{T}}y
\end{split}
$$

对所得到的估计进行分析，显然有：

$$
\begin{split}
E\hat{\beta}&=E\left[ \left( X_{c}^{\mathrm{T}}X_c \right) ^{-1}X_{c}^{\mathrm{T}}y \right] 
\newline
&=E\left[ \left( X_{c}^{\mathrm{T}}X_c \right) ^{-1}X_{c}^{\mathrm{T}}\left( 1\beta _0+X\beta +\epsilon \right) \right] 
\newline
&=\left( X_{c}^{\mathrm{T}}X_c \right) ^{-1}X_{}^{\mathrm{T}}\left( I-1\left( 1^{\mathrm{T}}1 \right) ^{-1}1^{\mathrm{T}} \right) \left( 1\beta _0+X\beta  \right) 
\newline
&=\left( X_{c}^{\mathrm{T}}X_c \right) ^{-1}X_{}^{\mathrm{T}}\left( I-1\left( 1^{\mathrm{T}}1 \right) ^{-1}1^{\mathrm{T}} \right) X\beta 
\newline
&=\left( X_{c}^{\mathrm{T}}X_c \right) ^{-1}X_{c}^{\mathrm{T}}X_c\beta 
\newline
&=\beta 
\end{split}
$$

这里利用了 $\left( I-1\left( 1^{\mathrm{T}}1 \right) ^{-1}1^{\mathrm{T}} \right) $ 对称幂等的性质。同理有

$$
\begin{split}
E\hat{\beta}_0&=E\left( \bar{y}-\bar{X}^{\mathrm{T}}\hat{\beta} \right) 
\newline
&=\left( 1^{\mathrm{T}}1 \right) ^{-1}1^{\mathrm{T}}Ey-\left( 1^{\mathrm{T}}1 \right) ^{-1}1^{\mathrm{T}}X\beta
\newline
&=\left( 1^{\mathrm{T}}1 \right) ^{-1}1^{\mathrm{T}}\left( 1\beta _0+X\beta \right) -\left( 1^{\mathrm{T}}1 \right) ^{-1}1^{\mathrm{T}}X\beta
\newline
&=\beta _0
\end{split}
$$

**因此 $\beta_0$，$\beta$ 都是无偏估计量。** 下面考虑他们的方差，则有

$$
\begin{split}
Var\left( \hat{\beta} \right) &=Var\left( \left( X_{c}^{\mathrm{T}}X_c \right) ^{-1}X_{c}^{\mathrm{T}}y \right) 
\newline
&=\left( X_{c}^{\mathrm{T}}X_c \right) ^{-1}X_{c}^{\mathrm{T}}Var\left( y \right) X_{c}^{}\left( X_{c}^{\mathrm{T}}X_c \right) ^{-1}
\newline
&=\sigma ^2\left( X_{c}^{\mathrm{T}}X_c \right) ^{-1}
\end{split}
$$

和

$$
\begin{split}
Var\left( \hat{\beta}_0 \right) &=Var\left( \bar{y}-\bar{X}^{\mathrm{T}}\hat{\beta} \right) 
\newline
&=Var\left[ \left( 1^{\mathrm{T}}1 \right) ^{-1}1^{\mathrm{T}}y-\left( 1^{\mathrm{T}}1 \right) ^{-1}1^{\mathrm{T}}X\left( X _{c}^{\mathrm{T}}X_c \right) ^{-1}X _{c}^{\mathrm{T}}y \right] 
\newline
&=\frac{\sigma ^2}{n^2}\left[ 1^{\mathrm{T}}1-0-0+1^{\mathrm{T}}X\left( X _{c}^{\mathrm{T}}X_c \right) ^{-1}X^{\mathrm{T}}1 \right]
\newline
&=\sigma ^2\left( \frac{1}{n}+\bar{X}^{\mathrm{T}}\left( X _{c}^{\mathrm{T}}X_c \right) ^{-1}\bar{X} \right) 
\end{split}
$$

上式只需注意到 $1^{\mathrm{T}}X_{c}=1^{\mathrm{T}}\left( I-1\left( 1^{\mathrm{T}}1 \right) ^{-1}1^{\mathrm{T}} \right) X=0$ 即可。最后，计算 $\beta_0$ 和 $\beta$ 之间的协方差：

$$
\begin{split}
Cov\left( \hat{\beta}_0,\hat{\beta} \right) &=Cov\left( \bar{y}-\bar{X}^{\mathrm{T}}\hat{\beta},\left( X _{c}^{\mathrm{T}}X_c \right) ^{-1}X _{c}^{\mathrm{T}}y \right)
\newline
&=\sigma ^2\left( 1^{\mathrm{T}}1 \right) ^{-1}1^{\mathrm{T}}\left( I-X\left( X _{c}^{\mathrm{T}}X_c \right) ^{-1}X _{c}^{\mathrm{T}} \right) X _{c}\left( X _{c}^{\mathrm{T}}X_c \right) ^{-1} 
\newline
&=\sigma ^2\left( 1^{\mathrm{T}}1 \right) ^{-1}1^{\mathrm{T}}\left( X _{c}\left( X _{c}^{\mathrm{T}}X _c \right) ^{-1}-X\left( X _{c}^{\mathrm{T}}X _c \right) ^{-1} \right) 
\newline
&=-\sigma ^2\left( 1^{\mathrm{T}}1 \right) ^{-1}1^{\mathrm{T}}X\left( X _{c}^{\mathrm{T}}X _c \right) ^{-1}
\newline
&=-\sigma ^2\bar{X}^{\mathrm{T}}\left( X _{c}^{\mathrm{T}}X _c \right) ^{-1}
\end{split}
$$


- **不带截距项的参数估计：**

假设 $\hat{\beta}$ 是参数 $\beta$ 的最小二乘估计量，$\hat{y}$ 是相应的 $y$ 的最小二乘估计，则根据最小二乘估计的思想，最小化以下目标函数即可：

$$
\hat{\beta}=\mathrm{arg}\ \mathop {\min} \limits_\beta\lVert y-X\beta \rVert _{2}^{2}
$$

令 $f\left( \beta \right) =\lVert y-X\beta \rVert _{2}^{2}$，对 $\beta$ 进行求导有：

$$
\begin{split}
\frac{\partial f\left( \beta \right)}{\partial \beta}&=-2X^{\mathrm{T}}y+2X^{\mathrm{T}}X\beta 
\newline
\frac{\partial ^2f\left( \beta \right)}{\partial \beta ^2}&=2X^{\mathrm{T}}X
\end{split}
$$

根据极值的必要条件，令上式为零有：

$$
\hat{\beta}=\left( X^{\mathrm{T}}X \right) ^{-1}X^{\mathrm{T}}y
$$

显然在求解中**并不需要** 假定 $\epsilon$ 一定服从正态分布。

对所得到的估计进行分析，显然有：

$$
\begin{split}
E\hat{\beta}&=E\left( X^{\mathrm{T}}X \right) ^{-1}X^{\mathrm{T}}y
\newline
&=\left( X^{\mathrm{T}}X \right) ^{-1}X^{\mathrm{T}}X\beta 
\newline
&=\beta 
\end{split}
$$

**因此 $\hat{\beta}$ 是无偏估计量。** 下面考虑他的方差：

$$
\begin{split}
Var\left( \hat{\beta} \right) &=Var\left( \left( X^{\mathrm{T}}X \right) ^{-1}X^{\mathrm{T}}y \right) 
\newline
&=\left( X^{\mathrm{T}}X \right) ^{-1}X^{\mathrm{T}}Var\left( y \right) X\left( X^{\mathrm{T}}X \right) ^{-1}
\newline
&=\sigma ^2\left( X^{\mathrm{T}}X \right) ^{-1}
\end{split}
$$


### 1.2 极大似然估计量
首先给出多元正态联合密度函数公式：

$$
f\left( x|\mu ,\Sigma \right) =\left[ \left( 2\pi \right) ^{\mathrm{T}}|\Sigma |^{-\frac{1}{2}} \right] \exp \left[ -\frac{1}{2}\left( x-\mu \right) ^{\mathrm{T}}\Sigma ^{-\frac{1}{2}}\left( x-\mu \right) \right]
$$

其中 $x$ 服从 $N(\mu, \Sigma)$ 分布。在本问题中，我们假定 $\Sigma =\sigma ^2I$，也就是随机变量是独立同分布于高斯分布的。

因此，根据极大似然原理得到似然函数：

$$
\mathscr{L}\left( \beta ,\sigma ^2 \right) =\left( 2\pi \sigma ^2 \right) ^{-\frac{n}{2}}\exp \left[ -\frac{1}{2\sigma ^2}\left( y-X\beta \right) ^{\mathrm{T}}\left( y-X\beta \right) \right]
$$

取对数得到对数似然函数：

$$
\ln \mathscr{L}\left( \beta ,\sigma ^2 \right) =-\frac{n}{2}\ln \left( 2\pi \sigma ^2 \right) -\frac{1}{2\sigma ^2}\left( y-X\beta \right) ^{\mathrm{T}}\left( y-X\beta \right) 
$$

上式分别对 $\beta$ 和 $\sigma^2$ 求偏导：

$$
\begin{split}
\frac{\partial \ln \mathscr{L}}{\partial \beta}&=-\frac{1}{2\sigma ^2}\left( -2X^{\mathrm{T}}y+X^{\mathrm{T}}X\beta \right) =0
\newline
\frac{\partial \ln \mathscr{L}}{\partial \sigma ^2}&=-\frac{n}{2\sigma ^2}+\frac{1}{2\sigma ^4}\left( y-X\beta \right) ^{\mathrm{T}}\left( y-X\beta \right) 
\end{split}
$$

显然根据第1式得到 $\beta$ 的极大似然估计量为：

$$
\hat{\beta}=\left( X^{\mathrm{T}}X \right) ^{-1}X^{\mathrm{T}}y
$$

并且我们还可以得到 $\sigma^2$ 的估计量：

$$
\hat{\sigma}^2=\frac{RSS}{n}
$$

其中 $RSS=\left( y-\hat{y} \right) ^{\mathrm{T}}\left( y-\hat{y} \right) $ 是残差平方和。

> 从这里可以看出，**在正态分布的假设下，最大似然法和最小二乘法的得到的估计量是一致的。** 因此，线性模型下使用最小二乘法或者极大似然法估计系数都可以。然而，这种等价性在[其他分布类型](/2018/06/glm/)的回归问题中不成立。


## 2 方差分析
### 2.1 带截距项的方差分解
首先给出模型总的离差平方和（$SYY$），回归平方和（$SS_{Reg}$）和残差平方和（$RSS$）的定义：

$$
\begin{split}
SYY&=\left( y-\bar{y}1 \right) ^{\mathrm{T}}\left( y-1\bar{y} \right) 
\newline
SS_{Reg}&=\left( \hat{y}-1\bar{y} \right) ^{\mathrm{T}}\left( \hat{y}-1\bar{y} \right) 
\newline
RSS&=\left( y-\hat{y} \right) ^{\mathrm{T}}\left( y-\hat{y} \right) 
\end{split}
$$

则对 $SYY$ 有：

$$
\begin{split}
SYY&=\left( y-\hat{y}+\hat{y}-1\bar{y} \right) ^{\mathrm{T}}\left( y-\hat{y}+\hat{y}-1\bar{y} \right) 
\newline
&=\left( y-\hat{y} \right) ^{\mathrm{T}}\left( y-\hat{y} \right) +\left( \hat{y}-1\bar{y} \right) ^{\mathrm{T}}\left( \hat{y}-1\bar{y} \right)+2\left( y-\hat{y} \right) ^{\mathrm{T}}\left( \hat{y}-1\bar{y} \right) 
\newline
&=SS_{\mathrm{Re}g}+RSS+2\left( \hat{y}-1\bar{y} \right) ^{\mathrm{T}}\left( y-\hat{y} \right) 
\end{split}
$$

考虑带有截距项的模型的参数估计：

$$
\begin{split}
\hat{\beta}_0&=\bar{y}-\bar{X}^{\mathrm{T}}\hat{\beta}
\newline
\hat{\beta}&=\left( X _{c}^{\mathrm{T}}X_c \right) ^{-1}X _{c}^{\mathrm{T}}y
\end{split}
$$

由此可以得到 $\hat{y}$ 为：

$$
\begin{split}
\hat{y}&=1\hat{\beta}_0+X\hat{\beta}
\newline
&=1\left[ \left( 1^{\mathrm{T}}1 \right) ^{-1}1^{\mathrm{T}}y-\left( 1^{\mathrm{T}}1 \right) ^{-1}1^{\mathrm{T}}X\left( X _{c}^{\mathrm{T}}X_c \right) ^{-1}X _{c}^{\mathrm{T}}y \right]+X\left( X _{c}^{\mathrm{T}}X_c \right) ^{-1}X _{c}^{\mathrm{T}}y
\newline
&=1\left( 1^{\mathrm{T}}1 \right) ^{-1}1^{\mathrm{T}}y+\left[ I-1\left( 1^{\mathrm{T}}1 \right) ^{-1}1^{\mathrm{T}} \right] X\left( X _{c}^{\mathrm{T}}X_c \right) ^{-1}X _{c}^{\mathrm{T}}y
\newline
&=1\left( 1^{\mathrm{T}}1 \right) ^{-1}1^{\mathrm{T}}y+X_c\left( X _{c}^{\mathrm{T}}X_c \right) ^{-1}X _{c}^{\mathrm{T}}y
\end{split}
$$

所以：

$$
\begin{split}
y-\hat{y}&=\left[ I-1\left( 1^{\mathrm{T}}1 \right) ^{-1}1^{\mathrm{T}} \right] y-X_c\left( X _{c}^{\mathrm{T}}X_c \right) ^{-1}X _{c}^{\mathrm{T}}y
\newline
\hat{y}-1\bar{y}&=X_c\left( X _{c}^{\mathrm{T}}X_c \right) ^{-1}X _{c}^{\mathrm{T}}y
\end{split}
$$

所以：

$$
\begin{split}
\left( \hat{y}-1\bar{y} \right) ^{\mathrm{T}}\left( y-\hat{y} \right) &=y^{\mathrm{T}}X_c\left( X _{c}^{\mathrm{T}}X_c \right) ^{-1}X _{c}^{\mathrm{T}}\left[ I-1\left( 1^{\mathrm{T}}1 \right) ^{-1}1^{\mathrm{T}} \right] y
\newline
&\quad-y^{\mathrm{T}}X_c\left( X _{c}^{\mathrm{T}}X_c \right) ^{-1}X _{c}^{\mathrm{T}}y
\end{split}
$$

而 $\left[ I-1\left( 1^{\mathrm{T}}1 \right) ^{-1}1^{\mathrm{T}} \right] $ 对称幂等，则有：

$$
\begin{split}
X _{c}^{\mathrm{T}}\left[ I-1\left( 1^{\mathrm{T}}1 \right) ^{-1}1^{\mathrm{T}} \right] &=X^{\mathrm{T}}\left[ I-1\left( 1^{\mathrm{T}}1 \right) ^{-1}1^{\mathrm{T}} \right] ^2
\newline
&=X^{\mathrm{T}}\left[ I-1\left( 1^{\mathrm{T}}1 \right) ^{-1}1^{\mathrm{T}} \right] 
\newline
&=X _{c}^{\mathrm{T}}
\end{split}
$$

所以得到 $\left( \hat{y}-1\bar{y} \right) ^{\mathrm{T}}\left( y-\hat{y} \right)$。这样我们就得到一个重要的**分解恒等式：**

$$
SYY=RSS+SS_{\mathrm{Reg}}
$$

**注意，在按照前面各个方差的定义下，只有带截距项的模型才使得上式成立。** 通常的模型都是带有截距的，因此在实际运用中，按照前面定义的不同离差是比较常见的。

> 在本节不同离差的定义下，不带截距项的模型分解恒等式不一定成立。


### 2.2 通用型方差分解式
除了以上形式的方差类型的定义方式外，还有一种模型方差分解的定义：

$$
\begin{split}
SYY&=y^{\mathrm{T}}y
\newline
SS_{Reg}&=\hat{y}^{\mathrm{T}}\hat{y}
\newline
RSS&=\left( y-\hat{y} \right) ^{\mathrm{T}}\left( y-\hat{y} \right) 
\end{split}
$$

将上式中的 $RSS$ 乘开得到：

$$
RSS=y^{\mathrm{T}}y+\hat{y}^{\mathrm{T}}\hat{y}-2y^{\mathrm{T}}\hat{y}
$$

考虑不带截距项的模型有 $\hat{\beta}=\left( X^{\mathrm{T}}X \right) ^{-1}X^{\mathrm{T}}y$，所以 $\hat{y}=X\left( X^{\mathrm{T}}X \right) ^{-1}X^{\mathrm{T}}y$，故而有：

$$
\begin{split}
\hat{y}^{\mathrm{T}}\hat{y}&=y^{\mathrm{T}}X\left( X^{\mathrm{T}}X \right) ^{-1}X^{\mathrm{T}}X\left( X^{\mathrm{T}}X \right) ^{-1}X^{\mathrm{T}}y
\newline
&=y^{\mathrm{T}}X\left( X^{\mathrm{T}}X \right) ^{-1}X^{\mathrm{T}}y
\newline
y^{\mathrm{T}}\hat{y}&=y^{\mathrm{T}}X\left( X^{\mathrm{T}}X \right) ^{-1}X^{\mathrm{T}}y=\hat{y}^{\mathrm{T}}\hat{y}
\end{split}
$$

所以将上式结果带入到 $RSS$ 中得到：

$$
\begin{split}
RSS&=y^{\mathrm{T}}y+\hat{y}^{\mathrm{T}}\hat{y}-2y^{\mathrm{T}}\hat{y}=y^{\mathrm{T}}y+\hat{y}^{\mathrm{T}}\hat{y}-2\hat{y}^{\mathrm{T}}\hat{y}
\newline
&=y^{\mathrm{T}}y-\hat{y}^{\mathrm{T}}\hat{y}
\newline
&=SYY-SS_{\mathrm{Re}g}
\end{split}
$$

这样就证明了方差分解公式成立。

> 在本节方差的定义下，带截距项模型的分解恒等式也成立。

通过上述分析可以看出，如果是带有截距项的模型，一般按照第一种定义进行分解；如果是不带截距项的模型，则按照式第二种定义进行分解。


### 2.3 方差分析表
- **带截距项的模型**

对带截距项的模型有：

$$
\begin{split}
\hat{\beta}_0&=\bar{y}-\bar{X}^{\mathrm{T}}\hat{\beta}
\newline
\hat{\beta}&=\left( X _{c}^{\mathrm{T}}X_c \right) ^{-1}X _{c}^{\mathrm{T}}y
\end{split}
$$

根据2.1节各个离差的定义，假定设计矩阵 $X$ 列满秩且 $\epsilon \sim N\left( 0,\sigma ^2I_n \right)$。

虑总离差平方和 $SYY$：

$$
\begin{split}
SYY&=\left( y-1\bar{y} \right) ^{\mathrm{T}}\left( y-1\bar{y} \right) 
\newline
&=y^{\mathrm{T}}\left( I-1\left( 1^{\mathrm{T}}1 \right) ^{-1}1^{\mathrm{T}} \right) ^{\mathrm{T}}\left( I-1\left( 1^{\mathrm{T}}1 \right) ^{-1}1^{\mathrm{T}} \right) y
\newline
&=y^{\mathrm{T}}\left( I-1\left( 1^{\mathrm{T}}1 \right) ^{-1}1^{\mathrm{T}} \right) y
\end{split}
$$

注意到 $P_1=1\left( 1^{\mathrm{T}}1 \right) ^{-1}1^{\mathrm{T}}$ 对称、幂等且 $I-P_1$ 秩为 $n-1$。对 $y=1\beta _0+X\beta+\epsilon $，容易得到 $y\sim N\left( \hat{\mu},\sigma ^2I_n \right) $，其中 $\hat{\mu}=1\beta _0+X\beta$，故而：

$$
\frac{SYY}{\sigma ^2}=\left( y/\sigma \right) ^{\mathrm{T}}\left( I-P_1 \right) \left( y/\sigma \right) \sim \chi _{n-1}^{2}\left( \frac{\hat{\mu}^{\mathrm{T}}\left( I-P_1 \right) \hat{\mu}}{\sigma ^2} \right) 
$$

考虑残差平方和 $RSS$：

$$
\begin{split}
RSS &=\left( y-\hat{y} \right) ^{\mathrm{T}}\left( y-\hat{y} \right) 
\newline
&=\left( y-1\hat{\beta}_0-X\hat{\beta} \right) ^{\mathrm{T}}\left( y-1\hat{\beta}_0-X\hat{\beta} \right) 
\newline
&=y^{\mathrm{T}}\left( I-P_1-P _{X_c}+P_1P _{X_c}+P _{X_c}P_1 \right) y
\end{split}
$$

其中 $P_{X_c}$ 秩为 $p$ 且满足：

$$
\begin{split}
P_{X_c}&=X_c\left( X _{c}^{\mathrm{T}}X_c \right) ^{-1}X _{c}^{\mathrm{T}}
\newline
&=\left( I-1\left( 1^{\mathrm{T}}1 \right) ^{-1}1^{\mathrm{T}} \right) X\left( X _{c}^{\mathrm{T}}X_c \right) ^{-1}X _{c}^{\mathrm{T}}
\newline
&=\left( I-P_1 \right) X\left( X _{c}^{\mathrm{T}}X_c \right) ^{-1}X _{c}^{\mathrm{T}}
\end{split}
$$

所以有 $\left( I-P_1 \right) P_{X_c}=P_{X_c}$，同理可以验证 $P_{X_c}\left( I-P_1 \right) =P_{X_c}$ 也成立。则对 $\left( I-P_1-P_{X_c}+P_1P_{X_c}+P_{X_c}P_1 \right) $ 有：

$$
\begin{split}
&I-P_1-P_{X_c}+P_1P_{X_c}+P_{X_c}P_1
\newline
&=I-\left( I-P_1 \right) P_{X_c}-P_1+P _{X_c}P_1
\newline
&=I-P_{X_c}-P_1+P_{X_c}P_1
\newline
&=I-P_1-P_{X_c}\left( I-P_1 \right) 
\newline
&=I-P_1-P_{X_c}
\end{split}
$$

可以验证 $I-P_1-P_{X_c}$ 是对称、幂等的且秩为 $n-p-1$，因此

$$
\frac{RSS}{\sigma ^2}=\left( \frac{y}{\sigma} \right) ^{\mathrm{T}}\left( I-P_1-P_{X_c} \right) \left( \frac{y}{\sigma} \right) \sim \chi _{n-p-1}^{2}
$$

这里的卡方分布中心没有偏移是因为

$$
\begin{split}
\hat{\mu}^{\mathrm{T}}\left( I-P_1-P_{X_c} \right) &=\left( 1\beta _0+X\beta \right) ^{\mathrm{T}}\left( I-P_1-P_{X_c} \right) 
\newline
&=\left( \beta _01^{\mathrm{T}}-\beta _01^{\mathrm{T}}P_1 \right) -\beta _01^{\mathrm{T}}P_{X_c}+\beta^{\mathrm{T}}X^{\mathrm{T}}-\beta^{\mathrm{T}}X^{\mathrm{T}}P_1-\beta^{\mathrm{T}}X^{\mathrm{T}}P_{X_c}
\newline
&=0-0+\beta^{\mathrm{T}}X^{\mathrm{T}}\left( I-P_1 \right) -\beta^{\mathrm{T}}X^{\mathrm{T}}P_{X_c}
\newline
&=0
\end{split}
$$

考虑回归平方和 $SS_{Reg}$：

$$
\begin{split}
SS_{Reg}&=\left( \hat{y}-1\bar{y} \right) ^{\mathrm{T}}\left( \hat{y}-1\bar{y} \right) 
\newline
&=y^{\mathrm{T}}X_c\left( X _{c}^{\mathrm{T}}X_c \right) ^{-1}X _{c}^{\mathrm{T}}y
\newline
&=y^{\mathrm{T}}P_{X_c}y
\end{split}
$$

结合 $y\sim N\left( \hat{\mu},\sigma ^2I_n \right)$，$\hat{\mu}=1\beta _0+X\beta+\mu$，故而

$$
\frac{SS_{Reg}}{\sigma ^2}=\left( y/\sigma \right) ^{\mathrm{T}}P_{X_c}\left( y/\sigma \right) \sim \chi _{p}^{2}\left( \frac{\tilde{\mu}^{\mathrm{T}}P _{X_c}\tilde{\mu}}{\sigma ^2} \right) 
$$

> 观察到 $\left( I-P_1-P_{X_c} \right) \times P_{X_c}=0$，故而 $RSS$ 和 $SS_{Reg}$ 是独立同分布的。

综合上述可得到带截距项的模型方差分析表

|方差来源 |平方和 |自由度 $df$ |均方 |
|:---|:---|:---|:---|
|回归 |$SS_{Reg}$ |$p$ |$SS_{Reg}/p$ |
|残差 |$RSS$ |$n-p-1$ |$RSS/(n-p-1)$ |
|总离差 |$SYY$ |$n-1$ | |

- **不带截距的模型**

对不带截距项的模型有：

$$
\hat{\beta}=\left( X^{\mathrm{T}}X \right) ^{-1}X^{\mathrm{T}}y
$$

根据式2.2各个离差的定义，假定设计矩阵 $X$ 列满秩且 $\epsilon \sim N\left( 0,\sigma ^2I_n \right) $。

考虑总离差平方和 $SYY$：

其中 $y=X\beta +\epsilon \sim \left( \tilde{\mu},\sigma ^2I_n \right) $，$\tilde{\mu}=X\beta $，故而

$$
\frac{SYY}{\sigma ^2}=\left( y/\sigma \right) ^{\mathrm{T}}\left( y/\sigma \right) \sim \chi _{n}^{2}\left( \frac{\tilde{\mu}^{\mathrm{T}}\tilde{\mu}}{\sigma ^2} \right) 
$$

考虑残差平方和 $RSS$

$$
\begin{split}
RSS=\left( y-\hat{y} \right) ^{\mathrm{T}}\left( y-\hat{y} \right) &=\left( y-X\hat{\beta} \right) ^{\mathrm{T}}\left( y-X\hat{\beta} \right) 
\newline
&=\left( y-X\left( X^{\mathrm{T}}X \right) ^{-1}X^{\mathrm{T}}y \right) ^{\mathrm{T}}\left( y-X\left( X^{\mathrm{T}}X \right) ^{-1}X^{\mathrm{T}}y \right) 
\newline
&=y^{\mathrm{T}}\left( I-P_X \right) y
\end{split}
$$

其中 $P_X=X\left( X^{\mathrm{T}}X \right) ^{-1}X^{\mathrm{T}}$ 是空间 $\mathscr{R}\left( X \right) $ 的正交投影算子。根据矩阵知识， $I-P_X$ 对称、幂等且秩为 $n-p$，所以

$$
\frac{RSS}{\sigma ^2}=\left( y/\sigma \right) ^{\mathrm{T}}\left( I-P_X \right) \left( y/\sigma \right) \sim \chi _{n-p}^{2}
$$

这里的卡方分布也是没有中心偏移的。

考虑回归平方和 $SS_{Reg}$

$$
SS_{Reg}=\hat{y}^{\mathrm{T}}\hat{y}=y^{\mathrm{T}}P_Xy
$$

结合 $y=X\beta +\epsilon \sim N\left( \tilde{\mu},\sigma ^2I_n \right)$，$\tilde{\mu}=X\beta +\mu $，故而

$$
\frac{SS_{Reg}}{\sigma ^2}=\left( y/\sigma \right) ^{\mathrm{T}}P_X\left( y/\sigma \right) \sim \chi _{p}^{2}\left( \frac{\tilde{\mu}^{\mathrm{T}}P_X\tilde{\mu}}{\sigma ^2} \right) 
$$

> 通过上述分析，观察到 $\left( I-P_X \right) \times P_X=0$，故而 $RSS$ 和 $SS_{Reg}$ 是独立同分布的。

最后，综合上述可得到不带截距项的模型方差分析表

|方差来源 |平方和 |自由度 $df$ |均方 |
|:---|:---|:---|:---|
|回归 |$SS_{Reg}$ |$p$ |$SS_{Reg}/p$ |
|残差 |$RSS$ |$n-p$ |$RSS/(n-p)$ |
|总离差 |$SYY$ |$n$ | |


## 3 统计显著性检验
### 3.1 模型显著性
- **带截距项的模型**

在这个检验前提下，我们的原假设和备择假设通常是

$$
\begin{split}
&H_0:\beta=0
\newline
&H_1:\beta\ne 0
\end{split}
$$

根据前述可知

$$
\begin{split}
&\frac{RSS}{\sigma^2}\ \sim \chi_{n-p-1}^{2}
\newline
&\frac{SS_{Reg}}{\sigma^2} \sim \chi_{p}^{2}\left( \frac{\hat{\mu}^{\mathrm{T}}P_{X_c}\hat{\mu}}{\sigma^2} \right) 
\end{split}
$$

其中 $\hat{\mu}=1\beta _0+X\beta$，故而当原假设成立时有 $\hat{\mu}=1\beta _0$，因此

$$
\begin{split}
\hat{\mu}^{\mathrm{T}}\left( I-P_1-P_{X_c} \right) \hat{\mu}&=\beta _{0}^{2}1^{\mathrm{T}}\left( I-P_1-P_{X_c} \right) 1
\newline
&=\beta _{0}^{2}\left( 1^{\mathrm{T}}1-1^{\mathrm{T}}1-1^{\mathrm{T}}P _{X_c}1 \right) 
\newline
&=-\beta _{0}^{2}1^{\mathrm{T}}P_{X_c}1
\end{split}
$$

而对 $1^{\mathrm{T}}P_{X_c}1$ 有

$$
\begin{split}
1^{\mathrm{T}}P_{X_c}1&=1^{\mathrm{T}}X_c\left( X_{c}^{\mathrm{T}}X_c \right) ^{-1}X_{c}^{\mathrm{T}}1
\newline
&=1^{\mathrm{T}}\left( I-1\left( 1^{\mathrm{T}}1 \right) ^{-1}1^{\mathrm{T}} \right) X\left( X_{c}^{\mathrm{T}}X_c \right) ^{-1}X_{c}^{\mathrm{T}}1
\newline
&=\left( 1^{\mathrm{T}}-1^{\mathrm{T}} \right) X\left( X_{c}^{\mathrm{T}}X_c \right) ^{-1}X_{c}^{\mathrm{T}}1
\newline
&=0
\end{split}
$$

所以 $\hat{\mu}^{\mathrm{T}}\left( I-P_1-P_{X_c} \right) \hat{\mu}=0$，同理可知 $\hat{\mu}^{\mathrm{T}}P_{X_c}\hat{\mu}=0$，故而在原假设成立的条件下有

$$
\begin{split}
\frac{RSS}{\sigma^2}\ &\sim \chi_{n-p-1}^{2}
\newline
\frac{SS_{Reg}}{\sigma^2} &\sim \chi_{p}^{2}
\end{split}
$$

所以，在 $H_0$ 成立的条件下有

$$
F_0=\frac{SS_{Reg}/p}{RSS/\left( n-p-1 \right)}\sim F_{p,n-p-1}
$$

这样当给定置信水平 $\alpha$，在原假设成立的条件下回归平方和 $SS_{Reg}$ 较小，而残差平方和 $RSS$则较大，因此检验统计量 $F_0$ 就应该较小。那么，拒绝原假设的条件就是 $F_0\ge F_{\alpha ,p,n-p-1}$。 这里 $F_{\alpha ,p,n-p-1}$ 表示上侧 $\alpha$ 分位数（下同）。


- **不带截距项的模型**

在这个检验前提下，我们的原假设和备择假设通常是

$$
\begin{split}
&H_0:\beta _1=\beta _2=\cdots =\beta _p=0
\newline
&H_1:\beta _j\ne 0,\mathrm{至少对一}个j\mathrm{成立}
\end{split}
$$

类似的，当原假设成立时有 $\hat{\mu}=0$，则

$$
\begin{split}
\frac{RSS}{\sigma^2}\ &\sim \chi_{n-p}^{2}
\newline
\frac{SS_{Reg}}{\sigma^2} &\sim \chi_{p}^{2}
\end{split}
$$

所以，在原假设成立的条件下有

$$
F_0=\frac{SS_{Reg}/p}{RSS/\left( n-p \right)}~F_{p,n-p}
$$

拒绝原假设的条件同前。

### 3.2 单个系数的检验

单个系数的检验通常和从模型中增删变量有关。一般来说，增加一个额外的变量到模型中去，不会使得回归平方和 减小，也不会使得残差平方和 $SS_{Reg}$ 增大。但是，回归平方和 $SS_{Reg}$ 的一点点增大能否充分的保证在模型中引入了一个冗余变量？我们之所以关心，是因为增加一个冗余变量到模型中确实会增大均方误差，因而将模型的可用性降低了。

检验任何一个单独的回归系数 $\beta_j$ 的假设通常是

$$
\begin{split}
&H_0:\beta _j=0
\newline
&H_1:\beta _j\ne 0
\end{split}
$$

- **带截距项的模型**

对带截距项的模型有

$$
\hat{\beta}=\left( X_{c}^{\mathrm{T}}X_c \right) ^{-1}X_{c}^{\mathrm{T}}y
$$

因为 $E\hat{\beta}=\beta$ 和 $Var\left( \hat{\beta} \right) =\sigma ^2\left( X_{c}^{\mathrm{T}}X_c \right) ^{-1}$，显然有

$$
\frac{\hat{\beta}_j-\beta _j}{\sqrt{\sigma ^2C _{jj}}}\sim N\left( 0,1 \right) 
$$

其中 $C_{jj}$ 是矩阵 $\left( X _{c}^{\mathrm{T}}X_c \right) ^{-1}$ 第 个对角元。但是由于一般情况下 $\sigma^2$ 未知，所以我们再考虑 $RSS$ 有

$$
\frac{RSS}{\sigma ^2}\ \sim \chi _{n-p-1}^{2}
$$

所以有

$$
\frac{\left( n-p-1 \right) s^2}{\sigma ^2}\sim \chi _{n-p-1}^{2}
$$

其中 $s^2$ 是样本的方差，是可以通过样本计算得到的。那么，在原假设成立的条件下有

$$
t_0 = \frac{\hat{\beta}_j}{\sqrt{s^2C _{jj}}} \sim t _{n-p-1}
$$

因此，给定置信水平 $\alpha$ 在原假设成立的条件下，上式 $t_0$ 的绝对值就不会特别大。那么，拒绝 $H_0$ 的条件就是：$| t_0 |>t_{\alpha /2,n-p-1}$。

对于截距项有

$$
\frac{\hat{\beta}_0-\beta _0}{\sqrt{\sigma ^2\left( 1/n+\bar{X}^{\mathrm{T}}\left( X _{c}^{\mathrm{T}}X_c \right) ^{-1}\bar{X} \right)}}\sim N\left( 0,1 \right) 
$$

因此

$$
t_0 = \frac{\hat{\beta}_j}{\sqrt{s^2(1/n + \bar{X} ^{\mathrm{T}}(X _c^{\mathrm{T}}X _c)^{-1}\bar{X})}} \sim t _{n-p-1}
$$

综述，给定置信水平 $\alpha$ 在原假设 $H_0$ 成立的条件下，拒绝 $H_0$ 的条件就是：$| t_0 |>t_{\alpha /2,n-p-1}$。

- **不带截距项的模型**

对不带截距项的模型有

$$
\hat{\beta}=\left( X_{}^{\mathrm{T}}X_{} \right) ^{-1}X_{}^{\mathrm{T}}y
$$

类似的，因为 $E\hat{\beta}=\beta $ 和 $E\hat{\beta}=\beta $，所以

$$
\frac{\hat{\beta}_j-\beta_j}{\sqrt{\sigma ^2C _{jj}}} \sim N\left( 0,1 \right) 
$$

其中 $C_{jj}$ 是矩阵 $(X^{\mathrm{T}}X)^{-1}$ 第 $j$ 个对角元。但是由于一般情况下 $\sigma^2$ 未知，所以我们再考虑 $RSS$ 有

$$
\frac{RSS}{\sigma ^2}\ \sim \chi _{n-p}^{2}
$$

所以有

$$
\frac{\left( n-p \right) s^2}{\sigma ^2} \sim \chi _{n-p}^{2}
$$

其中 $s^2$ 是样本的方差。那么在原假设成立的条件下有

$$
t_0 = \frac{\hat{\beta}_j}{\sqrt{s^2C _{jj}}} \sim t _{n-p}
$$

综述，给定置信水平 $\alpha$ 在原假设 $H_0$ 成立的条件下，拒绝 $H_0$ 的条件就是：$| t_0 |>t_{\alpha /2,n-p}$。


## 4 区间估计
### 4.1 单个回归系数的置信区间
- **带截距项的区间估计**

根据3.2节的分析和讨论，对于不带截距项的模型有

$$
\frac{\hat{\beta}_j-\beta _j}{\sqrt{s^2C _{jj}}}\sim t _{n-p-1}
$$

所以，给定置信水平 $\alpha$ 下 $\beta_j$ 置信区间为：

$$
\hat{\beta} _{j} - t _{\alpha /2,n-p-1}\sqrt{s^2C _{jj}}\le \beta _j\le \hat{\beta}_j+t _{\alpha /2,n-p-1}\sqrt{s^2C _{jj}}
$$


- **不带截距项的区间估计**

类似的，给定置信水平 $\alpha$ 下 $\beta_j$ 置信区间为

$$
\hat{\beta} _{j} - t _{\alpha /2,n-p}\sqrt{s^2C _{jj}}\le \beta _j\le \hat{\beta}_j+t _{\alpha /2,n-p}\sqrt{s^2C _{jj}}
$$

### 4.2 回归系数的联合置信区域
前面的置信区间都是针对单个系数进行的，因此置信水平只对一个区间有效。但是，很多问题中需要让置信水平对所有的区间有效，这就是联合置信区域（simultaneous confidence intervals）。

- **带截距项的模型**

由3.2节可知

$$
\left( \hat{\beta}-\beta \right) \sim N\left( 0,\sigma ^2\left( X_{c}^{\mathrm{T}}X_c \right) ^{-1} \right) 
$$

所以

$$
\frac{\left( \hat{\beta}-\beta \right) ^{\mathrm{T}}\left( X_{c}^{\mathrm{T}}X_c \right) ^{-1}\left( \hat{\beta}-\beta \right)}{\sigma ^2}\sim \chi _{p}^{2}
$$

这里 $\sigma^2$ 未知，为构造检验统计量再考虑 $RSS$ 有

$$
\frac{RSS}{\sigma ^2} \sim \chi _{n-p-1}^{2}
$$

所以得到

$$
\frac{\left( \hat{\beta}-\beta \right) ^{\mathrm{T}}\left( X_{c}^{\mathrm{T}}X_c \right) ^{-1}\left( \hat{\beta}-\beta \right) /p}{RSS/\left( n-p-1 \right)} \sim F_{p,n-p-1}
$$

那么在置信水平 $\alpha$ 下，所有的回归系数需要满足

$$
\frac{\left( \hat{\beta}-\beta \right) ^{\mathrm{T}}\left( X_{c}^{\mathrm{T}}X_c \right) ^{-1}\left( \hat{\beta}-\beta \right) /p}{RSS/\left( n-p-1 \right)}\le F_{\alpha ,p,n-p-1}
$$

- **不带截距项的模型**

同理可知，对于不带截距项的模型，在给定置信水平 $\alpha$ 下全部系数需要满足

$$
\frac{\left( \hat{\beta}-\beta \right) ^{\mathrm{T}}\left( X_{c}^{\mathrm{T}}X_c \right) ^{-1}\left( \hat{\beta}-\beta \right) /p}{RSS/\left( n-p-1 \right)}\le F_{\alpha ,p,n-p}
$$

对于前述联合置信区域的不等式描述的是一个**椭圆约束区域。** 当 $p=2$ 时，这个区域相当简单；但是当 $p\geq 2$ 时，问题就变得复杂的多。


## 5 模拟分析
我在[简单回归分析](/2019/07/simplelm/)中介绍了`R`求解线性模型的方法，这里我们根据前面的理论分析结果，用`R`编写系数求解的函数、单个系数的显著性检验函数。**这些代码附在最后供参考。**

我们考虑模型

$$
y=1\beta _0+X\beta+\epsilon 
$$

其中 $\beta _0=3$，$\beta=\left( 1,3,2,1,4,5,2,6 \right) ^{\mathrm{T}}$，$y\in \mathbb{R}^{100}$。根据前述分析，使用最小二乘法进行估计，在`R`软件中得到结果如下

![](/linearmodel/mylm-withb.png)

可以看出估计得结果很好，模型和单个系数都通过了检验。

当然，我们可以对相同的数据使用`R`软件求解带截距的线性模型，得到的结果如下图所示。可以看出，两者是一致的。

![](/linearmodel/rlm-withb.png)

此外，我们也可以对模型考虑不带截距项的估计。这样得到的结果如下图所示

![](/linearmodel/mylm-withoutb.png)

对比`R`软件内置函数求解的结果

![](/linearmodel/rlm-withoutb.png)

首先，两者得到的结果是一样的。**其次，我们发现对于一个本身带有截距项的模型使用不带截距项的方法求解，得到的效果就不好。** 因为可以看到变量 $x_4$ 的系数估计没有通过检验。

我们可以通过计算模型得到的残差，如学生化残差，通过分析残差的分布和对应的正态分布之间的异同，从而对回归结果做出分析和判断。如下图所示，是带截距模型的学生化残差分布图（左）以及不带截距的学生化残差分布图（右）

![](/linearmodel/studenterro.jpg)

从图上可以看出来，使用正确的模型得到的残差分布和正态分布较为接近，说明估计的效果好；而使用错误的模型得到的结果，残差分布和正态分布相差较大，估计的效果不好。


## 6 自编函数代码
### 6.1 模型求解
实际使用，主要运行该代码即可。这个代码整合了参数估计、检验等代码，将结果合并输出。

```r
##求解线性模型
mylm <- function(x,y,method = "ols",intercept = T){#
  # 本程序用来求解线性模型参数估计 区间估计和假设检验
  # method可选：ols,max 前者最小二乘 后者极大似然估计
  # intercept为T表示模型带截距 为F表示模型不带截距 

  # 求解估计量
  b <- mylmestimate(x,y,method = method,intercept = intercept)
  
  # 求解离差平方和分解
  rsquare <- mylmsquare(x,y,b,intercept = intercept)
  
  # 单个系数检验
  indivi <- mylmindivi(x,y,b,intercept = intercept)
  
  # 模型检验
  lmtest <- mylmtest(x,y,b,intercept = intercept)
  
  if(intercept){
    bname <- c("(Intercept)",paste0(rep("x",p),1:p))
  }
  else{
    bname <- paste0(rep("x",p),1:p)
  }
  
  estimate.test <- data.frame(
    "Beta" = bname,
    "Estimate" = b,
    "t.value" = indivi$Tt,
    "P.value" = indivi$Tp,
    "Sig." = indivi$Ts
  )
  
  # 输出结果
  if(intercept){
    cat("Coefficients:","","\n")
    print(estimate.test)
    cat("---","","\n")
    cat("Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1","","\n")
    cat("---","","\n")
    cat("Multiple R-squared:",(rsquare$sreg/rsquare$syy),"")
    cat("Adjusted R-squared:",(1-(rsquare$rss*(n-1))/(rsquare$syy*(n-p-1))),"\n")
    cat("F-statistic:",lmtest$f0,"")
    cat("on",p,"")
    cat("and",(n-p-1),"")
    cat("DF,","","")
    cat("p-value:",lmtest$p,"")
  }
  else{
    cat("Coefficients:","","\n")
    print(estimate.test)
    cat("---","","\n")
    cat("Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1","","\n")
    cat("---","","\n")
    cat("Multiple R-squared:",(rsquare$sreg/rsquare$syy),"")
    cat("Adjusted R-squared:",(1-(rsquare$rss*n)/(rsquare$syy*(n-p))),"\n")
    cat("F-statistic:",lmtest$f0,"")
    cat("on",p,"")
    cat("and",(n-p),"")
    cat("DF,","","")
    cat("p-value:",lmtest$p,"")
  }
}
```


### 6.2 参数估计

```r
##求解线性模型参数估计
mylmestimate <- function(x,y,method = "ols",intercept = T){#
  
  n <- length(y)
  
  if(method == "ols"){
    if(intercept){
      xbar0 <- apply(x,2,mean) # 按x的列计算均值
      xbar1 <- matrix(rep(xbar0,n),nrow = n,byrow = T)
      xcenter <- x - xbar1 # 得到按列中心化的矩阵x
      
      # 得到最小二乘估计
      b1 <- solve(t(xcenter)%*%xcenter)%*%t(xcenter)%*%y
      b0 <- mean(y) - t(xbar0)%*%b1
      
      b <- c(b0,b1)
    }
    else{
      # 得到最小二乘估计
      b <- solve(t(x)%*%x)%*%t(x)%*%y
    }
  }
  else{
    if(intercept){
      xbar0 <- apply(x,2,mean) # 按x的列计算均值
      xbar1 <- matrix(rep(xbar0,n),nrow = n,byrow = T)
      xcenter <- x - xbar1 # 得到按列中心化的矩阵x
      
      # 得到极大似然估计
      b1 <- solve(t(xcenter)%*%xcenter)%*%t(xcenter)%*%y
      b0 <- mean(y) - t(xbar0)%*%b1
      
      b <- c(b0,b1)
    }
    else{
      # 得到极大似然估计
      b <- solve(t(x)%*%x)%*%t(x)%*%y
    }
  }
  
  return(b)
}
```


### 6.3 假设检验

```r
##求解线性模型假设检验
mylmtest <- function(x,y,b,intercept = T){
  
  n <- dim(x)[1]
  p <- dim(x)[2]
  
  if(intercept){
    
    ybar <- mean(y)
    yhat <- b[1] + x%*%b[2:(p+1)]#计算回归值
    
    # 计算各个离差平方和
    rss <- t((y - yhat))%*%(y - yhat)
    sreg <- t((yhat - ybar))%*%(yhat - ybar)
    
    # 计算检验统计量
    f0 <- (sreg/p)/(rss/(n - p - 1))
    
    # 输出p值
    p <- 1 - pf(f0,p,n-p-1)
  }
  else{
    yhat <- x%*%b#计算回归值
    
    # 计算各个离差平方和
    rss <- t((y - yhat))%*%(y - yhat)
    sreg <- t(yhat)%*%(yhat)
    
    # 计算检验统计量
    f0 <- (sreg/p)/(rss/(n - p))
    
    # 输出p值
    p <- 1 - pf(f0,p,n-p)
  }
  
  modeltest <- data.frame("f0" = f0,"p" = p)
  
  return(modeltest)
}
```


### 6.4 单个系数检验

```r
##求解线性模型单个系数检验
mylmindivi <- function(x,y,b,intercept = T){
  
  n <- dim(x)[1]
  p <- dim(x)[2]
  
  if(intercept){
    
    yhat <- b[1] + x%*%b[2:(p+1)]#计算回归值
    s <- t((y - yhat))%*%(y - yhat)/(n-p-1) #计算MSE
    
    xbar0 <- apply(x,2,mean) # 按x的列计算均值
    xbar1 <- matrix(rep(xbar0,n),nrow = n,byrow = T)
    xcenter <- x - xbar1 # 得到按列中心化的矩阵x
    
    cc <- solve(t(xcenter)%*%xcenter)
    diagc <- diag(cc)
    
    Tt <- 1:(p+1);Tp <- Tt;Ts <- Tt
    Tt[1] <- b[1]/sqrt(s*(1/n + t(xbar0)%*%cc%*%xbar0))
    Tp[1] <- 2*(1 - pt(Tt[1],n-p-1))
    if(Tp[1] < 0.001 ){# 判断置信度
      Ts[1] <- c("***")
    }
    else if(Tp[1] < 0.01 && Tp[1] >= 0.001){
      Ts[1] <- c("**")
    }
    else if(Tp[1] < 0.05 && Tp[1] >= 0.01){
      Ts[1] <- c("*")
    }
    else if(Tp[1] < 0.1 && Tp[1] >= 0.05){
      Ts[1] <- c(".")
    }
    else{
      Ts[1] <- c(" ")
    }
    
    for(i in 2:(p+1)){
      Tt[i] <- b[i]/sqrt(s*diagc[i-1])
      Tp[i] <- 2*(1 - pt(Tt[i],n-p-1))
      if(Tp[i] < 0.001 ){# 判断置信度
        Ts[i] <- c("***")
      }
      else if(Tp[i] < 0.01 && Tp[i] >= 0.001){
        Ts[i] <- c("**")
      }
      else if(Tp[i] < 0.05 && Tp[i] >= 0.01){
        Ts[i] <- c("*")
      }
      else if(Tp[i] < 0.1 && Tp[i] >= 0.05){
        Ts[i] <- c(".")
      }
      else{
        Ts[i] <- c(" ")
      }
    }
  }
  else{
    yhat <- x%*%b#计算回归值
    s <- t((y - yhat))%*%(y - yhat)/(n-p) #计算MSE
      
    cc <- solve(t(x)%*%x)
    diagc <- diag(cc)
      
    Tt <- 1:p;Tp <- Tt;Ts <- Tt
    
    for(i in 1:p){
      Tt[i] <- b[i]/sqrt(s*diagc[i])
      Tp[i] <- 2*(1 - pt(Tt[i],n-p))
      if(Tp[i] < 0.001 ){# 判断置信度
        Ts[i] <- c("***")
        }
      else if(Tp[i] < 0.01 && Tp[i] >= 0.001){
        Ts[i] <- c("**")
      }
      else if(Tp[i] < 0.05 && Tp[i] >= 0.01){
        Ts[i] <- c("*")
      }
      else if(Tp[i] < 0.1 && Tp[i] >= 0.05){
        Ts[i] <- c(".")
      }
      else{
        Ts[i] <- c(" ")
      }
    }
  }
  
  indivitest <- data.frame("Tt" = Tt,"Tp" = Tp,"Ts" = Ts)
  return(indivitest)
}
```


### 6.5 区间估计

```r
##求解线性模型区间估计
mylminterval <- function(x,y,b,alpha = 0.05,intercept = T){
  
  n <- dim(x)[1]
  p <- dim(x)[2]
  
  if(intercept){
    
    Ti <- matrix(0,nrow = (p+1),ncol = 2)
    
    yhat <- b[1] + x%*%b[2:(p+1)]#计算回归值
    s <- t((y - yhat))%*%(y - yhat)/(n-p-1) #计算MSE
    
    xbar0 <- apply(x,2,mean) # 按x的列计算均值
    xbar1 <- matrix(rep(xbar0,n),nrow = n,byrow = T)
    xcenter <- x - xbar1 # 得到按列中心化的矩阵x
    
    cc <- solve(t(xcenter)%*%xcenter)
    diagc <- diag(cc)
    
    Ti[1,1] <- b[1] - qt((1-alpha/2),n-p-1)*sqrt(s*diagc[1])
    Ti[1,2] <- b[1] + qt((1-alpha/2),n-p-1)*sqrt(s*diagc[1])
    
    for(i in 2:(p+1)){
      
      Ti[i,1] <- b[i] - qt((1-alpha/2),n-p-1)*sqrt(s*diagc[1])
      Ti[i,2] <- b[i] + qt((1-alpha/2),n-p-1)*sqrt(s*diagc[1])
      
    }
  }
  else{
    
    Ti <- matrix(0,nrow = p,ncol = 2)
    
    yhat <- x%*%b#计算回归值
    s <- t((y - yhat))%*%(y - yhat)/(n-p) #计算MSE
    
    cc <- solve(t(x)%*%x)
    diagc <- diag(cc)
    
    for(i in 1:p){
      
      Ti[i,1] <- b[i] - qt((1-alpha/2),n-p-1)*sqrt(s*diagc[1])
      Ti[i,2] <- b[i] + qt((1-alpha/2),n-p-1)*sqrt(s*diagc[1])
      
    }
  }
  
  colnames(Ti) <- c("LowerBound","UpperBound")
  return(Ti)
}
```


### 6.6 离差平方和分解

```r
##求解线性模型离差平方和分解
mylmsquare <- function(x,y,b,intercept = T){
  if(intercept){
    
    ybar <- mean(y)
    yhat <- b[1] + x%*%b[2:(p+1)]#计算回归值
    
    # 计算各个离差平方和
    rss <- t((y - yhat))%*%(y - yhat)
    sreg <- t((yhat - ybar))%*%(yhat - ybar)
    syy <- t((y - ybar))%*%(y - ybar)
  }
  else{
    yhat <- x%*%b#计算回归值
    
    # 计算各个离差平方和
    rss <- t((y - yhat))%*%(y - yhat)
    sreg <- t(yhat)%*%(yhat)
    syy <- t(y)%*%y
  }
  
  rsquare <- data.frame("syy" = syy,"rss" = rss,"sreg" = sreg)
  return(rsquare)
}
```


### 6.7 各种残差

```r
##求解线性模型各种残差
mylmresidual <- function(x,y,b,intercept = T){
  
  n <- dim(x)[1]
  p <- dim(x)[2]
  
  if(intercept){
    
    yhat <- b[1] + x%*%b[2:(p+1)]# 计算回归值
    sig.hat <- t(y - yhat)%*%(y - yhat)/(n-p-1)# 计算方差估计值
     
    # 计算帽子矩阵对角元
    xbar0 <- apply(x,2,mean) # 按x的列计算均值
    xbar1 <- matrix(rep(xbar0,n),nrow = n,byrow = T)
    xcenter <- x - xbar1 # 得到按列中心化的矩阵x
    p1 <- matrix(1/n,nrow = n,ncol = n)
    h <- p1 + xcenter%*%solve(t(xcenter)%*%xcenter)%*%t(xcenter)
    diagh <- diag(h)
    
    # 计算原始残差
    originres <- y - yhat
    
    # 计算standardized residuals
    standres <- (y - yhat)/sqrt(sig.hat*rep(1,n))
    
    # 计算student residuals
    studres <- (y - yhat)/sqrt(sig.hat*(1 - diagh))
    
    # 计算PRESS residual
    pressres <- (y - yhat)/(1 - diagh)
  }
  else{
    yhat <- x%*%b# 计算回归值
    sig.hat <- t((y - yhat))%*%(y - yhat)/(n-p)# 计算方差估计值
    
    # 计算原始残差
    originres <- y - yhat
    
    # 计算帽子矩阵对角元
    h <- x%*%solve(t(x)%*%x)%*%t(x)
    diagh <- diag(h)
    
    # 计算standardized residuals
    standres <- (y - yhat)/sqrt(sig.hat*rep(1,n))
    
    # 计算student residuals
    studres <- (y - yhat)/sqrt(sig.hat*(1 - diagh))
    
    # 计算PRESS residual
    pressres <- (y - yhat)/(1 - diagh)
  }
  
  lmres <- data.frame("oringinal.Residual" = originres,
                      "Standardized.Residual" = standres,
                      "Student.Residual" = studres,
                      "PRESS.Residual" = pressres)
  return(lmres)
}
```