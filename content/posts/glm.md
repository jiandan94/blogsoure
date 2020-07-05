---
title: "广义线性模型及其一般求解方式"
date: 2018-06-15
toc: true
categories:
  - 笔记
tags:
  - 多元统计
  - 回归分析
---

## 1 指数族分布与广义线性模型
### 1.1 引入指数族分布
在线性模型中，一个重要的条件便是响应变量 $y$ 须服从正态分布。然而，实际问题中 $y$ 并不总是满足正态分布的假设。因此，我们考虑更加一般的指数族分布。

指数族分布定义如下：

$$
f\left( y;\theta ,\phi \right) =\exp \left[ \frac{y\theta -b\left( \theta \right)}{a\left( \phi \right)}+c\left( y,\phi \right) \right]
$$

其中 $a(\cdot)$，$b(\cdot)$ 和 $c(\cdot)$ 是已知给定的函数。参数 $\theta$ 为分布族的位置参数（location parameter），参数 $\phi$ 通常被称为分散参数（dispersion parameter）。一般来说，函数 $a(\phi)$ 通常有形式 $a(\phi)=\phi\cdot\omega$，其中 $\omega$ 是一个已知的常数。指数族分布包括常见的二项分布、泊松分布、正态分布和指数分布等。


### 1.2 联系回归问题
前面直接给出指数族分布，有点让人一时难以和回归问题构建联系。我们不妨回忆[线性模型](/2018/05/lm/)的内容，此时响应变量 $y$ 满足下面的正态分布

$$
N(x^T\beta, \sigma^2)
$$

我们写出它的密度函数具体表达式

$$
f\left( y \right) = \frac{1}{\sqrt{2\pi}\sigma}\exp\left[ -\frac{(y - x\beta)^2}{2\sigma^2} \right] = \exp\left[ -\frac{(y^2 - 2yx\beta + \beta^Tx^Tx\beta)}{2\sigma^2} - \frac{1}{2} \ln(2\pi\sigma^2) \right]
$$

整理一下就得到

$$
f\left( y \right) = \exp\left[ \frac{yx\beta - 1/2\beta^Tx^Tx\beta}{\sigma^2} - \frac{y^2}{2\sigma^2} - \frac{1}{2} \ln(2\pi\sigma^2) \right]
$$

**原来线性回归问题就是 $\theta = x\beta$ 和 $\phi = \sigma^2$ 的指数族问题。** 因此，我们保留线性模型中的线性结构假设，把正态分布约束推广成指数族分布，这就得到的更为一般的广义线性模型。

对应线性模型的假设，广义线性模型也有若干前提假设：

1.	观测 $y_1,\cdots ,y_n$ 是相互独立的，且对应的均值为 $\mu _1,\mu _2,\cdots ,\mu _n$；
2.	每个观测 $y_i$ 具有指数族分布；
3.	模型建立在线性预测因子 $\eta _1,\cdots ,\eta _n$ 上，其中 $\eta _i=x _{i}^{\mathrm{T}}\beta $，$x_i$ 是设计矩阵第 $i$ 个行向量；
4.	模型通过链接函数（link function）建立，其中：$\eta _i=g\left( \mu _i \right) ,\ i=1,2,\cdots ,n$；
5.	链接函数是单调可微的（其反函数存在）。

不难看出待求参数 $\beta$ 和指数族分布之间的关系链接：

$$
\beta \overset{\eta _i=g\left( \mu _i \right)}{\leftrightarrow}\mu _i\overset{\mu _i=\dot{b}\left( \theta _i \right)}{\leftrightarrow}\theta _i
$$

上述这种关系对我们逐步得到参数 $\beta$ 的求解公式意义重大。


## 2 广义线性模型的参数估计
在线性模型的假设下，[最小二乘法](/2019/08/forum-ols/)和极大似然法都能用于参数的求解。在广义线性模型中，我们无法写出二乘形式的优化函数。**因此，我们根据分布信息利用极大似然法来估计参数。** 模型的似然函数 $\mathscr{L}(\theta;Y)$ 为

$$
\mathscr{L}(\theta;Y) = \prod _{i=1}^n \exp \left[ \frac{y_i\theta_i -b\left( \theta-i \right)}{a\left( \phi_i \right)}+c\left( y_i,\phi_i \right) \right]
$$

为了方便理解算法的具体推导过程，下面首先介绍指数族的两个重要结论，然后再具体推导求解算法。


### 2.1 两个重要的结论
对于指数族分布，首先讨论两个重要的结论。根据密度函数可以得到相应的对数似然函数：

$$
\ln \mathscr{L}\left( \theta ;Y \right) =\sum_{i=1}^n{\frac{y_i\theta _i-b\left( \theta _i \right)}{a\left( \phi _i \right)}}+\sum_{i=1}^n{c\left( y_i,\phi _i \right)}
$$

这里仅假定 $y_i$ 是独立的。则似然函数对 $\theta$ 求导有：

$$
\frac{\partial \ln \mathscr{L}\left( \theta ;Y \right)}{\partial \theta}=\sum_{i=1}^n{\frac{y_i-\dot{b}\left( \theta _i \right)}{a\left( \phi _i \right)}}=\sum_{i=1}^n{S_i}
$$

似然函数对 $\theta$ 的二阶导为：

$$
\frac{\partial ^2\ln \mathscr{L}\left( \theta ;Y \right)}{\partial \theta ^2}=\sum_{i=1}^n{\frac{\ddot{b}\left( \theta _i \right)}{a\left( \phi _i \right)}}
$$

其中 $\dot{b}(\theta)$ 和 $\ddot{b}(\theta)$ 是 $b(\theta)$ 对 $\theta$ 的一阶导和二阶导（下同）。那么，大部分指数族分布（要求密度函数积分号和求导号可以交换顺序）满足：

$$
E\left( \frac{\partial \ln \mathscr{L}\left( \theta ;Y \right)}{\partial \theta} \right) =0, \quad E\left( \frac{\partial \ln \mathscr{L}\left( \theta ;Y \right)}{\partial \theta} \right) =0
$$

显然，根据上述两条性质可以推出：

$$
\mu _i=Ey_i=\dot{b}\left( \theta _i \right)
$$

和

$$
Var\left( y_i \right) =\ddot{b}\left( \theta _i \right) a\left( \phi \right) =\frac{\mathrm{d}\mu _i}{\mathrm{d}\theta _i}a\left( \phi \right) = Var _{\mu _i}a\left( \phi \right)
$$

### 2.2 极大似然估计
根据极大似然思想，对数似然函数对 $\beta$ 求导得到：

$$
\begin{split}
S\left( \beta \right) &=\frac{\partial \ln \mathscr{L}\left( \beta ;Y \right)}{\partial \beta}=\sum_{i=1}^n{\frac{\partial}{\partial \theta _i}\left( \frac{y_i\theta _i-b\left( \theta _i \right)}{a\left( \phi _i \right)} \right) \cdot \frac{\partial \theta _i}{\partial \mu _i}\cdot \frac{\partial \mu _i}{\partial \beta}}
\newline
&=\sum_{i=1}^n{\frac{y_i-\dot{b}\left( \theta _i \right)}{a\left( \phi _i \right)}\cdot \frac{\partial \theta _i}{\partial \mu _i}\cdot \frac{\partial \mu _i}{\partial \beta}}=\sum_{i=1}^n{\frac{y_i-\dot{b}\left( \theta _i \right)}{a\left( \phi _i \right)}\cdot \frac{1}{\frac{\partial \mu _i}{\partial \theta _i}}\cdot \frac{\partial \mu _i}{\partial \beta}}
\newline
&=\sum_{i=1}^n{\left( y_i-\dot{b}\left( \theta _i \right) \right) \frac{1}{a\left( \phi _i \right) \frac{\partial \mu _i}{\partial \theta _i}}\cdot \frac{\partial \mu _i}{\partial \beta}}=\sum_{i=1}^n{\left( y_i-\mu _i \right) \frac{1}{Var\left( y_i \right)}\cdot \frac{\partial \mu _i}{\partial \beta}}
\newline
&=\left( \frac{\partial \mu _1}{\partial \beta},\cdots ,\frac{\partial \mu _n}{\partial \beta} \right) \left( \begin{matrix}
	\frac{1}{Var\left( y_1 \right)}&		&		\newline
	&		\ddots&		\newline
	&		&		\frac{1}{Var\left( y_n \right)}\newline
\end{matrix} \right) \left( \begin{array}{c}
	y_1-\mu _1\newline
	\newline
	y_n-\mu _n\newline
\end{array} \right) 
\newline
&=\frac{\partial \mu ^{\mathrm{T}}}{\partial \beta}V^{-1}\left( y-\mu \right) =\left( \frac{\partial \mu}{\partial \beta ^{\mathrm{T}}} \right) ^{\mathrm{T}}V^{-1}\left( y-\mu \right) 
\newline
&=D^{\mathrm{T}}V^{-1}\left( y-\mu \right) 
\end{split}
$$

显然，这和之前求解得到的得分函数形式是一致的。为了得到具体的表达式，需要进一步求解 $D$ 的具体形式。

因为 $D=\frac{\partial \mu}{\partial \beta ^{\mathrm{T}}}$，而 $x_{i}^{\mathrm{T}}\beta =\eta _i=g\left( \mu _i \right) $，考虑到链接函数 $g(\cdot)$ 是单调可微函数，则其反函数存在，不妨设为 $h(\cdot)$，所以

$$
\frac{\partial \mu _i}{\partial \beta}=\frac{\partial}{\partial \beta}h\left( x _{i}^{\mathrm{T}}\beta \right) =\dot{h}\left( x _{i}^{\mathrm{T}}\beta \right) x_i
$$

进而

$$
\begin{split}
D&=\frac{\partial \mu}{\partial \beta ^{\mathrm{T}}}=\left( \frac{\partial \mu}{\partial \beta _1},\cdots ,\frac{\partial \mu}{\partial \beta _n} \right) =\left( \dot{h}\left( x _{1}^{\mathrm{T}}\beta \right) x_1,\cdots ,\dot{h}\left( x _{n}^{\mathrm{T}}\beta \right) x _n \right) 
\newline
&=\left( \begin{matrix}
	\dot{h}\left( x _{1}^{\mathrm{T}}\beta \right)&		&		\newline
	&		\ddots&		\newline
	&		&		\dot{h}\left( x _{n}^{\mathrm{T}}\beta \right)\newline
\end{matrix} \right) \left( x _1,\cdots ,x _n \right) 
\newline
&=\Delta X
\end{split}
$$

所以得到 $S(\beta)$ 的最终表达式为

$$
S\left( \beta \right) =\left( \Delta X \right) ^{\mathrm{T}}V^{-1}\left( y-\mu \right) =X^{\mathrm{T}}\Delta V^{-1}\left( y-\mu \right) 
$$

为了得到极大似然估计，我们接着求解对应的Fisher信息矩阵。根据定义有：

$$
I\left( \beta \right) =E\left( \frac{\partial}{\partial \beta ^{\mathrm{T}}}S\left( \beta \right) \right) =E\left( \frac{\partial ^2\ln \mathscr{L}}{\partial \beta ^{\mathrm{T}}\partial \beta} \right) 
$$

根据指数族的性质有

$$
E\left( \frac{\partial ^2\ln \mathscr{L}}{\partial \beta ^{\mathrm{T}}\partial \beta} \right) =-E\left( SS^{\mathrm{T}} \right) 
$$

考虑其第 $(i,j)$ 个元素 $S_{ij}=S_iS_j$，则

$$
\begin{split}
E\left( S_iS_j \right) &=E\left[ x _{i}^{\mathrm{T}}\Delta V^{-1}\left( y-\mu \right) x _{j}^{\mathrm{T}}\Delta V^{-1}\left( y-\mu \right) \right] 
\newline
&=E\left[ x _{i}^{\mathrm{T}}\Delta V^{-1}\left( y-\mu \right) \left( y-\mu \right) ^{\mathrm{T}}V^{-1}\Delta x _j \right] 
\newline
&=x _{i}^{\mathrm{T}}\Delta V^{-1}\Delta x _j
\end{split}
$$

所以

$$
I\left( \beta \right) =X^{\mathrm{T}}\Delta V^{-1}\Delta X
$$

这时，再回到头考虑得分函数，根据极值的必要条件有

$$
S\left( \hat{\beta} \right) =X^{\mathrm{T}}\Delta V^{-1}\left( y-\mu \right) =0
$$

将上式在真值 $\beta$ 处进行Taylor展开

$$
S\left( \hat{\beta} \right) =S\left( \beta \right) +\frac{\partial}{\partial \beta ^{\mathrm{T}}}S\left( \beta \right) \left( \hat{\beta}-\beta \right) =0
$$

所以

$$
\hat{\beta}=\beta +\left( -\frac{\partial S\left( \beta \right)}{\partial \beta ^{\mathrm{T}}} \right) ^{-1}S\left( \beta \right) 
$$

其中括号内的一项过于复杂，在实际求解中可以根据大数定律用其期望代替，由前面的讨论可知其期望即为 $I(\beta)$。这样，模型的极大似然估计可以近似成

$$
\begin{split}
\hat{\beta}&\doteq \beta +\left( X^{\mathrm{T}}\Delta V^{-1}\Delta X \right) ^{-1}X^{\mathrm{T}}\Delta V^{-1}\left( y-\mu \right) 
\newline
&=\beta +I\left( \beta \right) ^{-1}D^{\mathrm{T}}V^{-1}\left( y-\mu \right) 
\end{split}
$$


<!-- ### 2.3 渐近正态性
根据得到的估计量的近似表达式，其中等式右侧只有 $y$ 是服从正态分布的随机变量，那么得到的参数估计可以看成正态随机变量一个线性组合，显然估计量也服从正态分布，故而只需要求其均值和方差即可。

$\hat{\beta}$ 的均值

$$
E\hat{\beta}=E\left[ \beta +I\left( \beta \right) ^{-1}D^{\mathrm{T}}V^{-1}\left( y-\mu \right) \right] =\beta
$$

$\hat{\beta}$ 的方差

$$
\begin{split}
Var\left( \hat{\beta} \right) &=Var\left[ \beta +I\left( \beta \right) ^{-1}D^{\mathrm{T}}V^{-1}\left( y-\mu \right) \right] 
\newline
&=Var\left[ I\left( \beta \right) ^{-1}D^{\mathrm{T}}V^{-1}\left( y-\mu \right) \right] 
\newline
&=I\left( \beta \right) ^{-1}D^{\mathrm{T}}V^{-1}Var\left( y-\mu \right) V^{-1}DI\left( \beta \right) ^{-1}
\newline
&=I\left( \beta \right) ^{-1}\left( D^{\mathrm{T}}V^{-1}D \right) I\left( \beta \right) ^{-1}
\newline
&=I\left( \beta \right) ^{-1}
\end{split}
$$

所以，可知估计量服从 $N\left( \beta ,I\left( \beta \right) ^{-1} \right)$。有了渐近分布，我们就能构造显著性检验了。 -->


### 2.3 拟似然方法
广义线性模型中对于分布有前提假设，但是实际问题中并不知道具体的分布是什么。考虑到在求解问题中一个重要的信息就是信息函数，它和分布前两阶矩有关，故而可以假设分布的前两阶矩存在，从而推出不含分布信息的拟似然方法。

假定响应变量 $y_i$ 的均值为 $\mu_i$，方差函数为 $\mathrm{Var}y_i=a\left( \phi \right) \mathrm{Var} _{\mu _i}$，并且假定均值 $\mu_i$ 和 $x_i^T\beta$ 之间存在链结函数，链结函数的性质和广义线性模型的链结函数性质一样。据此，根据广义线性模型的思想实施拟似然方法。

构造逆得分函数

$$
S_i\left( \mu _i \right) =\frac{y_i-\mu _i}{\mathrm{Var}y_i}
$$

它满足

$$
\begin{cases}
	ES_i\left( \mu _i \right) =0&		\newline
	ES _{i}^{2}=E\left( -\frac{\partial S_i}{\partial \mu _i} \right)&		\newline
\end{cases}
$$

其中第二条性质是因为

$$
\begin{split}
E\left( -\frac{\partial S_i}{\partial \mu _i} \right) &=-E\left( \frac{\partial}{\partial \mu _i}\left( \frac{y_i-\mu _i}{\mathrm{Var}y_i} \right) \right) 
\newline
&=\frac{1}{\mathrm{Var}y_i}=E\left( S _{i}^{2} \right) 
\end{split}
$$

所以根据广义线性模型极大似然思想有

$$
\theta _i=\int _{y_i}^{\mu _i}{S_i\left( t \right) \mathrm{d}t}, \quad \theta \left( \mu \right) =\sum _{i=1}^n{\theta \left( \mu _i \right)}=\theta \left( \beta \right) 
$$

那么拟似然方法的得分函数可以写成

$$
S\left( \beta \right) =\frac{\partial}{\partial \beta}\theta \left( \beta \right) =\sum_{i=1}^n{\frac{\partial \theta _i}{\partial \beta}=\sum_{i=1}^n{\frac{\partial \mu _{i}^{\mathrm{T}}}{\partial \beta}\frac{\partial \theta _i}{\partial \mu _i}=D^{\mathrm{T}}V^{-1}\left( y-\mu \right)}}
$$

类似的可以得到Fisher信息矩阵为

$$
I\left( \beta \right) =D^{\mathrm{T}}V^{-1}D
$$


## 3 迭代求解算法
考虑模型的极大似然估计，其中括号内的一项过于复杂，在实际求解中可以根据大数定律用其期望代替，由前面的讨论可知其期望即为 $I(\beta)$。这样，模型的极大似然估计可以近似成

$$
\begin{split}
\hat{\beta}&\doteq \beta +\left( X^{\mathrm{T}}\Delta V^{-1}\Delta X \right) ^{-1}X^{\mathrm{T}}\Delta V^{-1}\left( y-\mu \right) 
\newline
&=\beta +I\left( \beta \right) ^{-1}D^{\mathrm{T}}V^{-1}\left( y-\mu \right)
\end{split} 
$$

这样就得到了估计量的求解迭代公式。


## 4 模拟分析
在`R`中，我们有`glm`函数求解广义线性模型。这里，我结合前面的分析自己编写的相应的参数估计和假设检验的函数，包括二项分布（逻辑回归）、Poisson分布。**这些代码附在最后供参考。**

以Poisson分布举例，考虑符合Poisson分布的观测，其中链接函数为

$$
\eta _i=\ln \mu _i
$$

我们生成 $\mathbb{R}^{1000\times 5}$ 的设计矩阵 $X$，系数 $\beta =\left( 1,3,2,4,5 \right) ^{\mathrm{T}}$，使用自编的广义线性模型求解算法求解该模型得到

![](/glm/myglm.png)

这里的`myglm`是我自编的广义线性模型求解函数，其中囊括了线性模型、二项分布模型（逻辑回归）和Poisson模型的参数估计和显著性检验。对比一下`R`中的`glm`包求解的结果

![](/glm/glm.png)

可以看到结果几乎是一样的。


## 5 自编函数代码
### 5.1 模型求解
实际使用，主要运行该代码即可。这个代码整合了参数估计、检验等代码，将结果合并输出。

```r
## This program is to solve the generalized linear models: logistic, poisson

myglm <- function(x,y,b0,alpha = 0.05,family = "poisson"){#
  
  options(digits = 4)
  
  # family has two choice: logistic and poisson
  # alpha is the significance of the interval estimate of the coefficients
  
  # solve the model
  if(family == "logistic"){
    glm.result <- myglmlogistic(x,y,beta1 = b0)
  }
  else if(family == "poisson"){
    glm.result <- myglmpoisson(x,y,beta1 = b0)
  }
  else{
    stop(cat("!!!!!!模型参数错误，请选择给定模型logistic,poisson之一"))
  }
  
  b <- glm.result$模型的解
  y.fit <- glm.result$拟合值
  num <- glm.result$算法迭代次数
  
  # hypothesis testing
  indivi <- myglmindivi(x,y,b,expr = family)
  
  # confidence interval
  inter.sig <- myglminterval(x,y,b,alpha = alpha,expr = family)
  
  bname <- paste0(rep("beta",length(b)),1:length(b))
  
  # the residuls
  # glm.residual <- y - y.fit
  # res1 <- summary(glm.residual)
  # names(res1) <- c("最小值","下四分位数","中位数","均值","上四分位数","最大值")
  
  estimate.test <- data.frame(
    "估计值" = b,
    "下界" = inter.sig[,2],
    "上界" = inter.sig[,3],
    "z值" = indivi$Tu,
    "P值(>|z|)" = indivi$Tp,
    "置信度" = indivi$Ts,check.names = F
  )
  
  # output the results
    cat("\n")
    cat("Call: 这是不带截距项的",family,"模型,","下面是模型的分析结果：",seq = "","\n")
    cat("\n")
    cat("参数估计结果:","\n")
    cat("\n")
    print(estimate.test)
    cat("---","","\n")
    cat("置信度:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1","","\n")
    cat("区间估计的置信水平为：",alpha,"\n")
    cat("---","\n")
    cat("Fisher信息矩阵的迭代次数为:",num,"次","\n")
}
```


### 5.2 Logistic模型的解

```r
##本程序用来求解广义线性模型中——logistic模型的解
myglmlogistic <- function(x,y,beta1,N = 5000,e = 1e-10){#
  # x is the design matrix
  # b is the starting value of the iteration
  # N is the upper bound of the times of the iteration
  # e is the convergence criteria
  
  n <- length(x)
  x <- as.matrix(x)
  y <- as.matrix(y)
  
  g <- expression(log(mu/(1 - mu)))# the link function
  h <- expression(1/(1 + exp(-eta)))# the inverse of the link funtion
  dh <- D(h,"eta")# the first derivative funtion of h
  b <- expression(log(1 + exp(theta)))# the b(theta) function of the pdf
  db <- D(b,"theta")# the first derivative funtion of b
  ddb <- D(db,"theta")# the second derivative funtion of b
  db.inverse <- expression(log(mu/(1 - mu)))# the inverse of db
  
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
    theta <- eval(db.inverse)
    v.inverse <- diag(1/as.vector(eval(ddb)))
    
    p1 <- solve(t(D)%*%v.inverse%*%D)
    p2 <- t(D)%*%v.inverse%*%delta
    p3 <- x%*%beta1 + v.inverse%*%(y - mu)
    
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
  eta <- x%*%beta1
  y.fit <- eval(h)# compute the fitting values of y
  glmlogistic.result <- list("模型的解" = beta1,"算法迭代次数" = k-1,"拟合值" = y.fit)
  return(glmlogistic.result)
}
```


### 5.3 Poisson模型的解

```r
##本程序用来求解广义线性模型中——Poisson模型的解
myglmpoisson <- function(x,y,beta1,N = 5000,e = 1e-10){#
	# x is the design matrix
	# b is the starting value of the iteration
	# N is the upper bound of the times of the iteration
	# e is the convergence criteria
	
	n <- length(x)
	x <- as.matrix(x)
	y <- as.matrix(y)
	
	g <- expression(log(mu))# the link function
	h <- expression(exp(eta))# the inverse of the link funtion
	dh <- D(h,"eta")# the first derivative funtion of h
	b <- expression(exp(theta))# the b(theta) function of the pdf
	db <- expression(exp(theta))# the first derivative funtion of b
	ddb <- expression(exp(theta))# the second derivative funtion of b
	db.inverse <- expression(log(mu))# the inverse of db
	
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
		theta <- eval(db.inverse)
		v.inverse <- diag(1/as.vector(eval(ddb)))
		
		p1 <- solve(t(D)%*%v.inverse%*%D)
		p2 <- t(D)%*%v.inverse%*%delta
		p3 <- x%*%beta1 + v.inverse%*%(y - mu)
		
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
	eta <- x%*%beta1
	y.fit <- eval(h)# compute the fitting values of y
	glmpoisson.result <- list("模型的解" = beta1,"算法迭代次数" = k-1,"拟合值" = y.fit)
	return(glmpoisson.result)
}
```


### 5.4 拟似然

```r
##
##This program is quasi-likelihood method
##

myquasimle <- function(x,y,beta1,family,N = 5000,e = 1e-10){#
  # x is the design matrix
  # b is the starting value of the iteration
  # N is the upper bound of the times of the iteration
  # e is the convergence criteria
  
  n <- length(x)
  x <- as.matrix(x)
  y <- as.matrix(y)
  
  v.mu <- expression(mu^2)# The covariance function of y
  
  if(family == "possion"){
    h <- expression(exp(eta))# the inverse of the link funtion
    dh <- D(h,"eta")# the first derivative funtion of h
  }
  else if(family == "logistic"){
    h <- expression(1/(1 + exp(-eta)))
    dh <- D(h,"eta")
  }
  else{
    print("分布参数错误, 请选择logistic或者possion分布！！！！")
  }
  
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
    v.inverse <- diag(1/as.vector(eval(v.mu)))
    
    p1 <- solve(t(D)%*%v.inverse%*%D)
    p2 <- t(D)%*%v.inverse%*%delta
    p3 <- x%*%beta1 + v.inverse%*%(y - mu)
    
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
  eta <- x%*%beta1
  y.fit <- eval(h)# compute the fitting values of y
  quasi.result <- list("模型的解" = beta1,"算法迭代次数" = k-1,"拟合值" = y.fit)
  return(quasi.result)
}
```


### 5.5 单个系数检验

```r
## 求解广义线性模型单个系数检验
myglmindivi <- function(x,y,b,expr){#
  # expr 表示模型：logistic, poisson
  # b 表示由广义线性模型估计得到的估计量
  n <- length(x)
  p <- length(b)
  
  beta.hat <- b
  D <- matrix(0, nrow = n, ncol = p)
  
  # get the D matrix
  if(expr == "logistic"){
    h <- expression(1/(1 + exp(-eta)))# the inverse of the link funtion
    dh <- D(h,"eta")# the first derivative funtion of h
    b <- expression(log(1 + exp(theta)))# the b(theta) function of the pdf
    db <- D(b,"theta")# the first derivative funtion of b
    ddb <- D(db,"theta")# the second derivative funtion of b
    db.inverse <- expression(log(mu/(1 - mu)))# the inverse of db
    
    eta <- x%*%beta.hat
    delta <- diag(as.vector(eval(dh)))
    D <- delta%*%x# the D matrix
    
    mu <- eval(h)
    theta <- eval(db.inverse)
    v.inverse <- diag(1/as.vector(eval(ddb)))# the inverse of the covariance matrix of y
  }
  else if(expr == "poisson"){
    h <- expression(exp(eta))# the inverse of the link funtion
    dh <- D(h,"eta")# the first derivative funtion of h
    db <- expression(exp(theta))# the first derivative funtion of b
    ddb <- expression(exp(theta))# the second derivative funtion of b
    db.inverse <- expression(log(mu))# the inverse of db
    
    eta <- x%*%beta.hat
    delta <- diag(as.vector(eval(dh)))
    D <- delta%*%x# the D matrix
    
    mu <- eval(h)
    theta <- eval(db.inverse)
    v.inverse <- diag(1/as.vector(eval(ddb)))# the inverse of the covariance matrix of y
  }
  else{
    stop(cat("!!!!!!模型参数错误，请选择给定模型logistic,poisson之一"))
  }
  
  cc <- solve(t(D)%*%v.inverse%*%D)
  diagc <- diag(cc)
  
  Tt <- 1:p;Tp <- Tt;Ts <- Tt
  
  for(i in 1:p){
    
    Tt[i] <- beta.hat[i]/sqrt(diagc[i])
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


### 5.6 区间估计

```r
## the confidence interval of the generalized linear model
myglminterval <- function(x,y,b,expr,alpha = 0.05){
  
  n <- length(x)
  p <- length(b)
  
  beta.hat <- b
  D <- matrix(0, nrow = n, ncol = p)
  
  # get the D matrix
  if(expr == "logistic"){
    h <- expression(1/(1 + exp(-eta)))# the inverse of the link funtion
    dh <- D(h,"eta")# the first derivative funtion of h
    b <- expression(log(1 + exp(theta)))# the b(theta) function of the pdf
    db <- D(b,"theta")# the first derivative funtion of b
    ddb <- D(db,"theta")# the second derivative funtion of b
    db.inverse <- expression(log(mu/(1 - mu)))# the inverse of db
    
    eta <- x%*%beta.hat
    delta <- diag(as.vector(eval(dh)))
    D <- delta%*%x# the D matrix
    
    mu <- eval(h)
    theta <- eval(db.inverse)
    v.inverse <- diag(1/as.vector(eval(ddb)))# the inverse of the covariance matrix of y
  }
  else if(expr == "poisson"){
    h <- expression(exp(eta))# the inverse of the link funtion
    dh <- D(h,"eta")# the first derivative funtion of h
    db <- expression(exp(theta))# the first derivative funtion of b
    ddb <- expression(exp(theta))# the second derivative funtion of b
    db.inverse <- expression(log(mu))# the inverse of db
    
    eta <- x%*%beta.hat
    delta <- diag(as.vector(eval(dh)))
    D <- delta%*%x# the D matrix
    
    mu <- eval(h)
    theta <- eval(db.inverse)
    v.inverse <- diag(1/as.vector(eval(ddb)))# the inverse of the covariance matrix of y
  }
  else{
    stop(cat("!!!!!!模型参数错误，请选择给定模型logistic, gompertz, weibull之一"))
  }
  
  Ti <- matrix(0,nrow = p,ncol = 2)
  
  
  cc <- solve(t(D)%*%v.inverse%*%D)
  diagc <- diag(cc)
  
  for(i in 1:p){
    Ti[i,1] <- beta.hat[i] - qnorm((1-alpha/2))*sqrt(diagc[1])
    Ti[i,2] <- beta.hat[i] + qnorm((1-alpha/2))*sqrt(diagc[1])  
  }
  
  out <- cbind(beta.hat,Ti)
  colnames(out) <- c("Estimator","LowerBound","UpperBound")
  return(out)
}
```