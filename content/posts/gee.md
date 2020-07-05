---
title: "广义估计方程简介"
date: 2018-06-22
toc: true
categories:
  - 笔记
tags:
  - 多元统计
  - 纵向数据分析
---

本文总结了广义估计方程中的参数估计，并推导了用于求解模型的Gauss—Newton迭代法。此外，针对每个内容，本文还给出了相应的`R`软件求解算法，并做了相应的模拟。


## 1 广义估计方程
假定观测数据中有 $y_1,\cdots ,y_n$ 个体，对于个体 $y_i$ 有 $y_{i1},\cdots ,y_{im}$ 个观测。其中对于 $y_1,\cdots ,y_n$ 个体之间是相互独立的，而每个个体内的观测 $y_{i1},\cdots ,y_{im}$ 则是具有相关性的。对于观测 $y_{ij}$ 对应有 $\eta_{ij}$，满足 $\eta _{ij}=x_{ij_1}\beta _1+\cdots +x _{ij_p}\beta _p$。相应的对 $y_{ij}$ 有假设 $\mu _{ij}=Ey _{ij}$，$\mathrm{Var}y_{ij}=a\left( \phi \right) \mathrm{Var}_{\mu _{ij}}$，并且 $\eta_{ij}$ 和 $\mu_{ij}$ 之间存在链接函数 $g(x)$ 使得 $\eta _{ij}=g\left( \mu _{ij} \right) =x _{ij}^{\mathrm{T}}\beta 
$。那么，估计参数 $\beta$ 就是广义估计方程中的一项重要任务。


## 2 IEE
由于纵向数据中组内观测之间存在相关性，所以不能直接离用[GLM](/2018/06/glm/)的结果。为了推导出GEE的参数估计的表达式，**首先假设纵向数据组内也是相互独立的。** 

> 这么做，是先通过一个特殊的清醒理解广义估计方程，然后再将得到的方法进行一般化推广，使之可以适用于广义估计方程。

另一方面，我们注意到在广义估计方程中没有分布的信息，只有均值和方差的假设。因而，极大似然法在这里不能直接使用。（没有线性关系的保证和正态分布假设，最小二乘显然也不能使用）

**这里，我们利用拟似然的方法来求解模型。** 似然方法只利用到了前二阶矩信息，那么在给定均值和方差的条件下，我们可以仿照极大似然的思路构造似然函数，这就是拟似然的思想。

首先构造 $S_{ij}\left( \mu _{ij} \right) $

$$
S_{ij}\left( \mu _{ij} \right) =\frac{y_{ij}-\mu_{ij}}{\mathrm{Var}y_{ij}}=\frac{y_{ij}-\mu _{ij}}{a\left( \phi \right) \mathrm{Var} _{\mu _{ij}}}
$$


那么相应的有

$$
\theta_{ij}\left( \mu_{ij} \right) =\int_{y_{ij}}^{\mu _{ij}}{S _{ij}\left( t \right) \mathrm{d}t}
$$

接着根据上述结果得到得分函数

$$
\begin{split}
S\left( \beta \right) &=\frac{\partial \theta \left( \mu \right)}{\partial \beta} = \sum_j{\sum_j \frac{\partial \theta_{ij}}{\partial \beta}} = \sum_j{\sum_j \frac{\partial \mu_{ij}}{\partial \beta}S_{ij}}
\newline
&= \sum_j{\sum_j \frac{\partial \mu_{ij}}{\partial \beta}\cdot \frac{y_{ij} - \mu_{ij}}{a(\phi)\mathrm{Var}\mu_{ij}}}
\newline
&= \sum_i {D_i^T V_i^{-1}(y_i - \mu_i)}
\end{split}
$$

其中 $D_{i}^{\mathrm{T}}=\Delta _iX_i$，这里 $X_i$ 表示第 $i$ 次的观测矩阵，$\Delta_i$ 则为对角阵

$$
\Delta _i=\left( \begin{matrix}
	\dot{h}\left( \eta _{i1} \right)&		&		&		\newline
	&		\dot{h}\left( \eta _{i2} \right)&		&		\newline
	&		&		\ddots&		\newline
	&		&		&		\dot{h}\left( \eta _{im} \right)\newline
\end{matrix} \right) 
$$


## 3 GEE
根据上述求解的思想，对于GEE我们有

$$
S\left( \beta \right) =\sum_i{D_{i}^{\mathrm{T}}\left( \mathrm{Cov}y_i \right) ^{-1}\left( y_i-\mu _i \right)}
$$

其中 $\mathrm{Cov}y_i=A_{i}^{1/2}R_iA_{i}^{1/2}$。显然，$\beta$ 要计算参数 我们需要知道 $R_i$（组内相关结构），显然问题中 $R_i$ 未知，故而为了估计得以进行须给定 $R_{i0}$，称之为working correlation matrix。

这样，我们令

$$
D=\left( \begin{array}{c}
	D_1\newline
	\vdots\newline
	D_n\newline
\end{array} \right) ,
\quad
V=\left( \begin{matrix}
	\mathrm{Cov}y_1&		&		\newline
	&		\ddots&		\newline
	&		&		\mathrm{Cov}y_n\newline
\end{matrix} \right),
\quad
y=\left( \begin{array}{c}
	y_1\newline
	\vdots\newline
	y_n\newline
\end{array} \right),
\quad
\mu =\left( \begin{array}{c}
	\mu _1\newline
	\vdots\newline
	\mu _n\newline
\end{array} \right) 
$$

那么 $S(\beta)$ 可以写成

$$
S\left( \beta \right) =D^{\mathrm{T}}V^{-1}\left( y-\mu \right) 
$$

则相应的加权迭代最小二乘法可以写成

$$
\beta =\left( D^{\mathrm{T}}V^{-1}D \right) ^{-1}D^{\mathrm{T}}V^{-1}\Delta \left[ \Delta ^{-1}\left( y-\mu \right) +X\beta \right] 
$$


## 4 数据模拟
根据前面的分析可以编写求解模型的`R`代码，**这些都附在最后。**

本实验中我们考虑10个观测个体，每个个体观测100次，共有1000次观测。数据按照如下方式产生：设计矩阵 $X$ 是服从 $[0,0.1]$ 均匀分布，系数给定为 $\beta =\left( 1,3,2,4,5 \right) ^{\mathrm{T}}$，链接函数 $g(x) = \ln(x)$，组内协方差矩阵先生成 $100\times100$ 的服从 $[0,1]$ 均匀分布的矩阵，在将之和其转置相乘，取乘积的0.4倍作为组内协方差矩阵。这样，使用GEE模型，初始迭代向量取 $\left( 1.2,3.5,2.1,4.2,5.5 \right)$ 进行求解模型得到如下结果

![](/gee/gee.png)

可见最终的效果还是可以的。


## 5 程序代码

```r
##This program is quasi-likelihood method
mygee <- function(x,y,beta1,v.inverse,N = 5000,e = 1e-10){#
  # x is the design matrix
  # b is the starting value of the iteration
  # N is the upper bound of the times of the iteration
  # e is the convergence criteria
  # v.inverse is the inverse of the covariance function of y
  
  n <- length(x)
  x <- as.matrix(x)
  y <- as.matrix(y)
  
  h <- expression(exp(eta))# the inverse of the link funtion
  dh <- D(h,"eta")# the first derivative funtion of h
  
  k <- 1
  beta0 <- beta1 + 1
  while(sum((beta0 - beta1)^2) >= e){
    
    beta0 <- beta1
    
    # compute the D matrix
    eta <- x%*%beta1
    delta <- diag(as.vector(eval(dh)))
    D <- delta%*%x
    
    # compute v.inverse—the inverse matrix of the covariance matrix of y
    mu <- eval(h)
    
    p1 <- solve(t(D)%*%v.inverse%*%D)
    p2 <- t(D)%*%v.inverse%*%delta
    p3 <- x%*%beta1 + solve(delta)%*%(y - mu)
    
    beta1 <- p1%*%p2%*%p3# the kth estimates
    
    # check if it is divergent
    if(k > N){
      cat("算法不收敛，已达到最大迭代次数：",N,"\n")
      cat("此时得到的解为：","","\n")
      print(beta1)
      break
    }
    else{
      k <- k + 1
    }
  }
  
  colnames(beta1) <- c("估计值")
  rownames(beta1) <- paste0(rep("系数",length(beta1)),1:length(beta1))
  # eta <- x%*%beta1
  # y.fit <- eval(h)# compute the fitting values of y
  quasi.result <- list("模型的解" = beta1,"算法迭代次数" = k-1)
  return(quasi.result)
}
```