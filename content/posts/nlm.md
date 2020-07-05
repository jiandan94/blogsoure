---
title: "非线性模型的参数估计和统计性质"
date: 2018-05-18
toc: true
categories:
  - 笔记
tags:
  - 多元统计
  - 回归分析
---

这个笔记总结了非线性模型的极大似然估计量和渐进性质，并推导了用于求解模型的Gauss—Newton迭代法。此外，针对每个内容，我还给出了相应的`R`软件求解算法，并做了相应的模拟。


## 1 非线性模型
[线性模型](/2018/05/lm/)建立在自变量和响应变量之间呈线性关系的基础上，但实际数据并不总是如此。当我们没有额外信息认为两者之间的关系为线性时，非线性模型便成了一种选择。

考虑非线性模型

$$
y=f\left( X,\beta \right) +\epsilon 
$$

其中 $y\in \mathbb{R}^n$，$X\in \mathbb{R}^{n\times p}$，$n$ 表示观测数，$p$ 表示变量数。$\beta \in \mathbb{R}^n$ 表示回归系数。$\epsilon$ 一般假设成为独立同分布的 $N\left( 0,\sigma ^2I_n \right) $ 高斯随机变量。


## 2 模型的极大似然估计
如果我们假设模型的随机误差项是独立同分于均值为0方差为常数 $\sigma^2$ 的正态分布，那么可以考虑极大似然估计法来估计模型的解。

首先，根据模型和假设得到似然函数

$$
\mathscr{L}\left( \beta ,\sigma ^2 \right) =\frac{1}{\left( 2\pi \sigma ^2 \right) ^{n/2}}\exp \left[ -\frac{1}{2\sigma ^2}\left( y-f\left( X,\beta \right) \right) ^{\mathrm{T}}\left( y-f\left( X,\beta \right) \right) \right] 
$$

对上式取对数，得到对数似然函数

$$
\ln \mathscr{L}\left( \beta ,\sigma ^2 \right) =-\frac{n}{2}\ln \left( 2\pi \sigma ^2 \right) -\frac{1}{2\sigma ^2}\left( y-f\left( X,\beta \right) \right) ^{\mathrm{T}}\left( y-f\left( X,\beta \right) \right) 
$$

将对数似然方程对 $\beta$ 进行求偏导，根据极值必要条件令为零，从而得到得分方程（score equation）

$$
\begin{split}
\frac{\partial \ln \mathscr{L}\left( \beta ,\sigma ^2 \right)}{\partial \beta}&=\frac{\partial}{\partial \beta}\left[ -\frac{1}{2\sigma ^2}\left( y-f\left( X,\beta \right) \right) ^{\mathrm{T}}\left( y-f\left( X,\beta \right) \right) \right] 
\newline
&=-\frac{1}{2\sigma ^2}\frac{\partial \left( y-f\left( X,\beta \right) \right) ^{\mathrm{T}}}{\partial \beta}\left( y-f\left( X,\beta \right) \right) 
\newline
&=-\frac{1}{\sigma ^2}\frac{\partial f\left( X,\beta \right) ^{\mathrm{T}}}{\partial \beta}\left( y-f\left( X,\beta \right) \right) 
\newline
&=-\frac{1}{\sigma ^2}\left( \frac{\partial f\left( X,\beta \right)}{\partial \beta ^{\mathrm{T}}} \right) ^{\mathrm{T}}\left( y-f\left( X,\beta \right) \right) 
\newline
&=0
\end{split}
$$

上式利用了矩阵微分中的结论

$$
\frac{\mathrm{d}f\left( Y\left( X \right) \right)}{\mathrm{d}X}=\frac{\mathrm{d}Y^{\mathrm{T}}}{\mathrm{d}X}\cdot \frac{\mathrm{d}f\left( Y \right)}{\mathrm{d}Y}
$$

其中 $X\in \mathbb{R}^{n\times 1}$，$Y\in \mathbb{R}^{n\times 1}$，$f\in \mathbb{R}$ 的结论。进一步，上式等价于

$$
\left( \frac{\partial f\left( X,\beta \right)}{\partial \beta ^{\mathrm{T}}} \right) ^{\mathrm{T}}\left( y-f\left( X,\beta \right) \right) =0
$$

由于 $f(X, \beta)$ 表示均值函数，将其记作 $\mu$，并用 $\hat{\mu}$ 表示 $f(X,\beta)$ 中 $\beta$ 被其估计量 $b$ 替换的结果。那么根据上式可以推导出模型的估计量满足

$$
\left( \frac{\partial \hat{\mu}}{\partial \beta ^{\mathrm{T}}} \right) ^{\mathrm{T}}\left( y-\hat{\mu} \right) =0
$$

这时我们记 $\frac{\partial \hat{\mu}}{\partial \beta ^{\mathrm{T}}}=D\in \mathbb{R}^{n\times p}$，其中 $D_{ij}=\frac{\partial f\left( X_i,\beta \right)}{\partial \beta _j}$，那么就有

$$
D^{\mathrm{T}}\left( y-\hat{\mu} \right) =0
$$

> **这里可以看出极大似然估计量和最小二乘估计量的形式是一致的。** 之所以使用似然法，是因为在线性模型、[广义线性模型](/2018/06/glm/)中似然法更普适。

## 3 迭代求解算法
在非线性回归模型中，求解最小二乘估计量（极大似然估计量）的一个广泛应用的方法是将期望函数线性化，然后利用Gauss—Newton迭代法进行求解。

考虑非线性模型，将其在 点处进行Taylor展开

$$
y=f\left( X,\beta \right) +\epsilon =f\left( X,b_0 \right) +\frac{\partial f\left( X,\beta \right)}{\partial \beta ^{\mathrm{T}}}\left( \beta -b_0 \right) +\epsilon 
$$

令 $y_0=y-f\left( X,b_0 \right) =y-f_0$，$\theta _0=\beta -b_0$，$D_0=\left[ \frac{\partial f\left( X,\beta \right)}{\partial \beta ^{\mathrm{T}}} \right] _{\beta =b_0}$，则上式可以写成

$$
y_0=D_0\theta _0+\epsilon 
$$

那么 $\theta_0$ 的最小二乘估计量为

$$
\begin{split}
\hat{\theta}_0 &=\left( D _{0}^{\mathrm{T}}D_0 \right) ^{-1}D _{0}^{\mathrm{T}}y_0
\newline
&=\left( D _{0}^{\mathrm{T}}D_0 \right) ^{-1}D _{0}^{\mathrm{T}}\left( y-f_0 \right)
\end{split} 
$$

因为 $\theta _0=\beta -b_0$，我们用

$$
b_1=b_0+\hat{\theta}_0
$$

作为未知参数 $\beta$ 的一个修正估计。那么，我们就按照这个逻辑不停得修正 $\beta$，也就是说有

$$
\begin{split}
b_{k+1}&=b_k+\hat{\theta}_k
\newline
&=b_k+\left( D_{k}^{\mathrm{T}}D_k \right) ^{-1}D_{k}^{\mathrm{T}}\left( y-f_k \right) 
\end{split}
$$

当这种修正直到收敛（前后两个估计的改变量非常小）时，即

$$
\frac{\lVert b_{k+1}-b_k \rVert}{\lVert b_k \rVert}<\delta 
$$

迭代结束，其中 $\delta$ 是某个很小的数，比如说 $10^{-6}$。


## 4 估计量的渐近正态性
根据前面分析可知 $\beta$ 估计量为

$$
\hat{\beta}\doteq \beta +\left( D^{\mathrm{T}}D \right) ^{-1}D^{\mathrm{T}}\left( y-f \right) 
$$

其中 $D=\frac{\partial \mu}{\partial \beta ^{\mathrm{T}}}=\frac{\partial f\left( X,\beta \right)}{\partial \beta ^{\mathrm{T}}}$。对上式进行分析，等式右侧的随机项只有 $y$ 一项，这时我们考察有

$$
\begin{split}
E\hat{\beta}&=\beta +E\left[ \left( D^{\mathrm{T}}D \right) ^{-1}D^{\mathrm{T}}\left( y-f \right) \right] 
\newline
&=\beta +\left( D^{\mathrm{T}}D \right) ^{-1}D^{\mathrm{T}}E\left( y-f \right) 
\newline
&=\beta +\left( D^{\mathrm{T}}D \right) ^{-1}D^{\mathrm{T}}\left( Ey-f \right) 
\newline
&=\beta
\end{split} 
$$

和

$$
\begin{split}
Var\left( \hat{\beta} \right) &=Var\left[ \beta +\left( D^{\mathrm{T}}D \right) ^{-1}D^{\mathrm{T}}\left( y-f \right) \right] 
\newline
&=Var\left[ \left( D^{\mathrm{T}}D \right) ^{-1}D^{\mathrm{T}}\left( y-f \right) \right] 
\newline
&=\left( D^{\mathrm{T}}D \right) ^{-1}D^{\mathrm{T}}Var\left( y-f \right) D\left( D^{\mathrm{T}}D \right) ^{-1}
\newline
&=\sigma ^2\left( D^{\mathrm{T}}D \right) ^{-1}
\end{split}
$$

进一步，上式求出了得分函数

$$
S\left( \beta \right) =\frac{1}{\sigma ^2}D^{\mathrm{T}}\left( y-\mu \right) 
$$

这是一个 $\mathbb{R}^{p\times 1}$ 中的列向量，并且

$$
ES=\frac{1}{\sigma ^2}ED^{\mathrm{T}}\left( y-\mu \right) =D^{\mathrm{T}}\left( Ey-\mu \right) =0
$$

所以

$$
\begin{split}
Cov\left( S \right) &=\frac{1}{\sigma ^4}E\left( S-ES \right) \left( S-ES \right) ^{\mathrm{T}}=\frac{1}{\sigma ^4}ESS^{\mathrm{T}}
\newline
&=\frac{1}{\sigma ^4}E\left[ D^{\mathrm{T}}\left( y-\mu \right) \left( y-\mu \right) ^{\mathrm{T}}D \right] 
\newline
&=\frac{1}{\sigma ^4}D^{\mathrm{T}}E\left[ \left( y-\mu \right) \left( y-\mu \right) ^{\mathrm{T}} \right] D
\newline
&=\frac{1}{\sigma ^4}D^{\mathrm{T}}Var\left( y \right) D=\frac{1}{\sigma ^2}D^{\mathrm{T}}D
\end{split}
$$

上式得到的协方差矩阵称之为Fisher信息矩阵，记做 $I(\beta)$。所以，$\hat{\beta}$ 的协方差矩阵也可以写成

$$
Cov\left( \beta \right) =I^{-1}\left( \beta \right) 
$$

重新考察 $\hat{\beta}$ 的近似表达式，由于 $\hat{\beta}$ 是近似关于 $y$ 的一个线性组合，而 $y$ 是服从正态分布的，那么 $\hat{\beta}$ 也近似服从正态分布，根据所求的结果有

$$
\hat{\beta}\ \sim N\left( \beta ,I^{-1}\left( \beta \right) \right) 
$$

**如果 $y$ 不是独立的，协方差矩阵为 $V$，可以推出下面结果**

$$
S\left( \beta \right) =D^{\mathrm{T}}V^{-1}\left( y-\mu \right) 
$$

且

$$
\hat{\beta}\ \dot{\sim} N\left( \beta ,I^{-1}\left( \beta \right) \right) 
$$


## 5 模拟和实例
根据前面的分析可以自己用`R`编写非线性模型求解函数，以及显著性检验的函数。**相关代码附在最后供参考。**

> 迭代算法毕竟是局部最优，模拟的初始值最好取在真实值附近。**如果是实际问题，可以根据经验、利用OLS估计等最为初始值。**

### 5.1 Logistic模型

考虑Logistic增长模型

$$
y=\frac{\beta _1}{1+\beta _2e^{-\beta _3x}}+\varepsilon 
$$

在本例中取 $\beta_1=5$，$\beta_2=6$，$\beta_3 = 3$；并且 服从标准正态分布，样本观测数取5000。 $\epsilon$ 是服从标准正态分布的噪声。在模拟中，使用Gauss—Newton迭代法给定初始迭代值 $b_0=\left( 4.5,5.2,3.5 \right) ^{\mathrm{T}}$，程序迭代5次后收敛，得到系数估计如下图所示

![](/nlm/nlm-logistic.png)

可见算法的效果还是不错的。


### 5.2 Gompertz模型
考虑Gompertz模型

$$
y=\beta _1\exp \left( -\beta _2e^{-\beta _3x} \right) +\varepsilon 
$$

在本例中取 $\beta_1=30$，$\beta_2=14$，$\beta_3 = 3$；并且 $x$ 服从标准正态分布，样本观测数取5000。 $\epsilon$ 是服从标准正态分布的噪声。在模拟中，使用Gauss—Newton迭代法给定初始迭代值 $b_0=\left( 24.5,15.2,3.5 \right) ^{\mathrm{T}}$，程序迭代6次后收敛，得到系数估计如下图所示

![](/nlm/nlm-gompertz.png)

可见算法的效果还是不错的。


### 5.3 Weibull模型
考虑Weibull模型

$$
y=\beta _1-\beta _2\exp \left( -\beta _3x^{\beta _4} \right) +\varepsilon 
$$

在本例中取 $\beta_1=5$，$\beta_2=1$，$\beta_3 = 3$，$\beta_4=2$；并且 $x$ 服从 $[0,1]$ 的均匀分布，样本观测数取5000。$\epsilon$ 是服从标准正态分布的噪声。在模拟中，使用Gauss—Newton迭代法给定初始迭代值 ，程序迭代6次后收敛，得到系数估计如下图所示

![](/nlm/nlm-weibull.png)

可见算法的效果还是不错的。


### 5.4 Michaelis-Menten模型
考虑Michaelis-Menten模型

$$
y=\frac{\beta _1x}{\beta _2+x}+\varepsilon 
$$

对puromycin数据

|Concentration ($y$)|	Velocity ($x_1$)|Velocity ($x_2$) |
|:--- |:--- |:--- |
|0.02|47|76|
|0.06|97|107|
|0.11|123|139|
|0.22|152|159|
|0.56|191|201|
|1.10|200|207|

考虑上述Michaelis-Menten模型，使用Gauss—Newton迭代法给定初始迭代值 $b_0=\left( 205.00,0.08 \right) ^{\mathrm{T}}$，程序迭代5次后收敛，得到系数估计

$$
b=\left( 212.6837,0.0641 \right) ^{\mathrm{T}}
$$

利用得到的参数估计计算拟合值，画出图像如下所示

![](/nlm/nlm-puromycin.png)

图中圆点表示真实数据，虚线是得到的Michaelis-Menten模型估计方程。从图中看出，拟合的效果还是不错的。


### 5.5 渐近正态性验证
考虑5.1中的Logistic模型，在本次试验中分别取观测数目为50,500,1000和5000来进行模拟，每次模拟都做500次，然后画出每次参数估计中 $\beta_1$ 的核密度曲线，并与对应的正态密度曲线进行对比，结果如下所示

![](/nlm/nlm-anormal.png)

从图中可以看出，$\beta_1$ 的渐近正态性得到了验证。


## 6 自编函数代码
### 6.1 Logistic模型系数求解

```r
mynlmlogistic <- function(x,y,b,N = 5000,e = 1e-10){#
  # x is the design matrix
  # b is the starting value of the iteration
  # N is the upper bound of the times of the iteration
  # e is the convergence criteria
  
  n <- length(x)

  mylogistic <- expression(a1/(1 + a2*exp(-a3*t)))
  partial.a1 <- D(mylogistic,"a1")
  partial.a2 <- D(mylogistic,"a2")
  partial.a3 <- D(mylogistic,"a3")
  
  k <- 1
  bk <- b
  b <- b + 1# to enter the iteration
  while(sum((bk - b) ^ 2) / sum(b ^ 2) >= e) {
    # get the D matrix
    D <- matrix(0, nrow = n, ncol = 3)
    a1 <- bk[1]
    a2 <- bk[2]
    a3 <- bk[3]
    for (i in 1:n) {
      t <- x[i]
      D[i, 1] <- eval(partial.a1)
      D[i, 2] <- eval(partial.a2)
      D[i, 3] <- eval(partial.a3)
    }
    
    t <- x
    f <- eval(mylogistic)# compute the value of f
    b <- bk# denote the (k)th valuve by b
    bk <- bk + solve(t(D)%*%D)%*%t(D)%*%(as.matrix(y - f))# denote the (k+1)th value by bk
    
    # check if it is divergent
    if(k > N){
      cat("算法不收敛，已达到最大迭代次数：",N,"\n")
      cat("此时得到的解为：","","\n")
      print(bk)
      break
    }
    else{
      k <- k + 1
    }
  }
  
  colnames(bk) <- c("估计值")
  rownames(bk) <- paste0(rep("系数",length(bk)),1:length(bk))
  logistic.result <- list("模型的解" = bk,"算法迭代次数" = k-1)
  return(logistic.result)
}
```


### 6.2 Gompertz模型系数求解

```r
mynlmgompertz <- function(x,y,b,N = 5000,e = 1e-10){#
  # x is the design matrix
  # b is the starting value of the iteration
  # N is the upper bound of the times of the iteration
  # e is the convergence criteria
  
  n <- length(x)
 
  mygompertz <- expression(a1*exp(-a2*exp(-a3*t)))
  partial.a1 <- D(mygompertz,"a1")
  partial.a2 <- D(mygompertz,"a2")
  partial.a3 <- D(mygompertz,"a3")
  
  k <- 1
  bk <- b
  b <- b + 1# to enter the iteration
  while(sum((bk - b) ^ 2) / sum(b ^ 2) >= e) {
    # get the D matrix
    D <- matrix(0, nrow = n, ncol = 3)
    a1 <- bk[1]
    a2 <- bk[2]
    a3 <- bk[3]
    for (i in 1:n) {
      t <- x[i]
      D[i, 1] <- eval(partial.a1)
      D[i, 2] <- eval(partial.a2)
      D[i, 3] <- eval(partial.a3)
    }
    t <- x
    f <- eval(mygompertz)# compute the value of f
    b <- bk# denote the (k)th valuve by b
    bk <- bk + solve(t(D)%*%D)%*%t(D)%*%(as.matrix(y - f))# denote the (k+1)th value by bk
    
    # check if it is divergent
    if(k > N){
      cat("算法不收敛，已达到最大迭代次数：",N,"\n")
      cat("此时得到的解为：","","\n")
      print(bk)
      break
    }
    else{
      k <- k + 1
    }
  }
  
  colnames(bk) <- c("估计值")
  rownames(bk) <- paste0(rep("系数",length(bk)),1:length(bk))
  gompertz.result <- list("模型的解" = bk,"算法迭代次数" = k-1)
  return(gompertz.result)
}
```


### 6.3 Weibull模型系数求解

```r
mynlmweibull <- function(x,y,b,N = 5000,e = 1e-10){#
  # x is the design matrix
  # b is the starting value of the iteration
  # N is the upper bound of the times of the iteration
  # e is the convergence criteria
  
  n <- length(x)
  
  myweibull <- expression(a1 - a2*exp(-a3*(t^a4)))
  partial.a1 <- D(myweibull,"a1")
  partial.a2 <- D(myweibull,"a2")
  partial.a3 <- D(myweibull,"a3")
  partial.a4 <- D(myweibull,"a4")
  
  k <- 1
  bk <- b
  b <- b + 1# to enter the iteration
  while(sum((bk - b) ^ 2) / sum(b ^ 2) >= e) {
    # get the D matrix
    D <- matrix(0, nrow = n, ncol = 4)
    a1 <- bk[1]
    a2 <- bk[2]
    a3 <- bk[3]
    a4 <- bk[4]
    for (i in 1:n) {
      t <- x[i]
      D[i, 1] <- eval(partial.a1)
      D[i, 2] <- eval(partial.a2)
      D[i, 3] <- eval(partial.a3)
      D[i, 4] <- eval(partial.a4)
    }
    
    t <- x
    f <- eval(myweibull)# compute the value of f
    b <- bk# denote the (k)th valuve by b
    bk <- bk + solve(t(D)%*%D)%*%t(D)%*%(as.matrix(y - f))# denote the (k+1)th value by bk
    
    # check if it is divergent
    if(k > N){
      cat("算法不收敛，已达到最大迭代次数：",N,"\n")
      cat("此时得到的解为：","","\n")
      print(bk)
      break
    }
    else{
      k <- k + 1
    }
  }
  
  colnames(bk) <- c("估计值")
  rownames(bk) <- paste0(rep("系数",length(bk)),1:length(bk))
  weibull.result <- list("模型的解" = bk,"算法迭代次数" = k-1)
  return(weibull.result)
}
```


### 6.4 Michaelis-Menten模型系数求解

```r
mynlmmicmen <- function(x,y,b,N = 5000,e = 1e-10){#
  # x is the design matrix
  # b is the starting value of the iteration
  # N is the upper bound of the times of the iteration
  # e is the convergence criteria
  
  n <- length(x)
  
  mymicmen <- expression(a1*t/(a2 + t))
  partial.a1 <- D(mymicmen,"a1")
  partial.a2 <- D(mymicmen,"a2")
  
  k <- 1
  bk <- b
  b <- b + 1# to enter the iteration
  while(sum((bk - b) ^ 2) / sum(b ^ 2) >= e) {
    # get the D matrix
    D <- matrix(0, nrow = n, ncol = 2)
    a1 <- bk[1]
    a2 <- bk[2]
    for (i in 1:n) {
      t <- x[i]
      D[i, 1] <- eval(partial.a1)
      D[i, 2] <- eval(partial.a2)
    }
    
    t <- x
    f <- eval(mymicmen)# compute the value of f
    b <- bk# denote the (k)th valuve by b
    bk <- bk + solve(t(D)%*%D)%*%t(D)%*%(as.matrix(y - f))# denote the (k+1)th value by bk
    
    # check if it is divergent
    if(k > N){
      cat("算法不收敛，已达到最大迭代次数：",N,"\n")
      cat("此时得到的解为：","","\n")
      print(bk)
      break
    }
    else{
      k <- k + 1
    }
  }
  
  colnames(bk) <- c("估计值")
  rownames(bk) <- paste0(rep("系数",length(bk)),1:length(bk))
  mymicmen.result <- list("模型的解" = bk,"算法迭代次数" = k-1)
  return(mymicmen.result)
}
```


### 6.5 单个系数检验

```r
## 求解非线性模型单个系数检验
mynlmindivi <- function(x,y,b,expr){#
  # expr 表示模型：logistic, gompertz, weibull
  n <- length(x)
  p <- length(b)
  
  t <- x
  D <- matrix(0, nrow = n, ncol = p)
  
  # get the D matrix
  if(expr == "logistic"){
    mylogistic <- expression(a1/(1 + a2*exp(-a3*t)))
    partial.a1 <- D(mylogistic,"a1")
    partial.a2 <- D(mylogistic,"a2")
    partial.a3 <- D(mylogistic,"a3")
    a1 <- b[1]
    a2 <- b[2]
    a3 <- b[3]
    D[, 1] <- eval(partial.a1)
    D[, 2] <- eval(partial.a2)
    D[, 3] <- eval(partial.a3)
  }
  else if(expr == "gompertz"){
    mygompertz <- expression(a1*exp(-a2*exp(-a3*t)))
    partial.a1 <- D(mygompertz,"a1")
    partial.a2 <- D(mygompertz,"a2")
    partial.a3 <- D(mygompertz,"a3")
    a1 <- b[1]
    a2 <- b[2]
    a3 <- b[3]
    D[, 1] <- eval(partial.a1)
    D[, 2] <- eval(partial.a2)
    D[, 3] <- eval(partial.a3)
  }
  else if(expr == "weibull"){
    myweibull <- expression(a1 - a2*exp(-a3*(t^a4)))
    partial.a1 <- D(myweibull,"a1")
    partial.a2 <- D(myweibull,"a2")
    partial.a3 <- D(myweibull,"a3")
    partial.a4 <- D(myweibull,"a4")
    a1 <- b[1]
    a2 <- b[2]
    a3 <- b[3]
    a4 <- b[4]
    D[, 1] <- eval(partial.a1)
    D[, 2] <- eval(partial.a2)
    D[, 3] <- eval(partial.a3)
    D[, 4] <- eval(partial.a4)
  }
  else{
    stop(cat("!!!!!!模型参数错误，请选择给定模型logistic, gompertz, weibull之一"))
  }
  
  cc <- solve(t(D)%*%D)
  diagc <- diag(cc)
  
  Tt <- 1:p;Tp <- Tt;Ts <- Tt
  
  for(i in 1:p){
    
    Tt[i] <- b[i]/sqrt(diagc[i])
    Tp[i] <- 2*(1 - pnorm(Tt[i]))
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
  
  indivitest <- data.frame("Tu" = Tt,"Tp" = Tp,"Ts" = Ts)
  return(indivitest)
}
```


### 6.6 区间估计

```r
##求解非线性模型区间估计
mynlminterval <- function(x,y,b,expr,alpha = 0.05){
  
  n <- length(x)
  p <- length(b)
  
  t <- x
  D <- matrix(0, nrow = n, ncol = p)
  
  # get the D matrix
  if(expr == "logistic"){
    mylogistic <- expression(a1/(1 + a2*exp(-a3*t)))
    partial.a1 <- D(mylogistic,"a1")
    partial.a2 <- D(mylogistic,"a2")
    partial.a3 <- D(mylogistic,"a3")
    a1 <- b[1]
    a2 <- b[2]
    a3 <- b[3]
    D[, 1] <- eval(partial.a1)
    D[, 2] <- eval(partial.a2)
    D[, 3] <- eval(partial.a3)
  }
  else if(expr == "gompertz"){
    mygompertz <- expression(a1*exp(-a2*exp(-a3*t)))
    partial.a1 <- D(mygompertz,"a1")
    partial.a2 <- D(mygompertz,"a2")
    partial.a3 <- D(mygompertz,"a3")
    a1 <- b[1]
    a2 <- b[2]
    a3 <- b[3]
    D[, 1] <- eval(partial.a1)
    D[, 2] <- eval(partial.a2)
    D[, 3] <- eval(partial.a3)
  }
  else if(expr == "weibull"){
    myweibull <- expression(a1 - a2*exp(-a3*(t^a4)))
    partial.a1 <- D(myweibull,"a1")
    partial.a2 <- D(myweibull,"a2")
    partial.a3 <- D(myweibull,"a3")
    partial.a4 <- D(myweibull,"a4")
    a1 <- b[1]
    a2 <- b[2]
    a3 <- b[3]
    a4 <- b[4]
    D[, 1] <- eval(partial.a1)
    D[, 2] <- eval(partial.a2)
    D[, 3] <- eval(partial.a3)
    D[, 4] <- eval(partial.a4)
  }
  else{
    stop(cat("!!!!!!模型参数错误，请选择给定模型logistic, gompertz, weibull之一"))
  }
  
  Ti <- matrix(0,nrow = p,ncol = 2)
  
  
  cc <- solve(t(D)%*%D)
  diagc <- diag(cc)
  
  for(i in 1:p){
    Ti[i,1] <- b[i] - qnorm((1-alpha/2))*sqrt(diagc[1])
    Ti[i,2] <- b[i] + qnorm((1-alpha/2))*sqrt(diagc[1])  
  }
  
  out <- cbind(b,Ti)
  colnames(out) <- c("Estimator","LowerBound","UpperBound")
  return(out)
}
```