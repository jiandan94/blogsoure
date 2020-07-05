---
title: "提升R的运算速度"
date: 2020-03-23
toc: true
categories:
  - R语言
tags:
  - R语言
  - 算法
---


最近有读者私信希望介绍**R并行计算**的方法。在处理较大规模问题时，默认单线程的R显得力不从心。提升运算速度确实在数据分析中意义重大，其中酸爽**调参侠**定深有感触。

不过我觉得饭要一口一口吃，如果算法先天不足、设计不合理，指望通过某个函数或程序包暴力提升运算速度，就很难得到满意的结果。

因此，我打算从一些简单的细节入手，给大家分享R运算速度提升的常见方法和思维误区。

## 1. 内置函数与向量化思维
R对常见的运算如矩阵乘积，进行过专门的优化。那么，我们应当有意识将算法中的相关部分转用R优化过的内置函数处理，从而降低计算成本。

又如，我们在学习模型时，常常看到作者用分量形式展示模型的求解过程——虽然这令人易懂，但让我们在设计算法时极易忽视将其向量化，最后造成计算开销骤增。

### - 案例展示

我们不妨考虑不带截距项的单变量最小二乘模型

$$\hat{\beta}=arg\min\sum_{i=1}^n\left( y_i - x_i\hat{\beta}\right)^2$$

很容易得到最小二乘解为

$$\hat{\beta} = \frac{\sum_{i=1}^nx_iy_i}{\sum_{i=1}^nx_i^2}$$

这时，你可能将程序写成

```r{.line-numbers}
## 算法1
ols1 <- function(x, y){
  a1 <- 0
  a2 <- 0
  for(i in 1:length(y)){
    a1 <- a1 + x[i]*y[i]
    a2 <- a2 + x[i]^2
  }
  return(a1/a2)
}
```

假如你知道R中`*`的用法，你可能会将程序写成

```r{.line-numbers}
## 算法2
ols2 <- function(x, y){
  return(sum(x*y)/sum(x^2))
}
```

**当然你还可以写成向量形式**

```r{.line-numbers}
## 算法3
ols3 <- function(x, y){
  return(c(x%*%y/x%*%x))
}
```

为了揭示这三者计算效率的差别，我们不妨从标准正态分布中产生1000个观测，重复计算100次的结果如下

```r
> x <- rnorm(1000)
> y <- x*2 + 0.1*rnorm(1000)
> system.time(# 算法1
+   replicate(100, ols1(x, y))
+ )
用户 系统 流逝 
0.02 0.00 0.02 
> system.time(# 算法2
+   replicate(100, ols2(x, y))
+ )
用户 系统 流逝 
0.01 0.00 0.01 
> system.time(# 算法3
+   replicate(100, ols3(x, y))
+ )
用户 系统 流逝 
   0    0    0 
```

结果是显然的：向量化的运算比`for`循环效率明显高。因此，我们要记住 **能避免使用`for`循环就避免使用，** 因为`for`循环在R中很低效。

善用R内置函数与向量化思维是很多新手易忽视的技巧。这个例子提醒大家须将其铭记于心。


## 2. 警惕apply族函数

稍有R常识的读者可能怀疑我题目写错了，因为大家印象中`apply`族函数可以简化代码和提升效率，所以要广泛使用而不是警惕啊。简化代码我同意，但提升效率真的不一定正确。

### - apply函数

我们先看最常见的`apply`函数的用法

```r
apply(X, MARGIN, FUN, ...)
```

其中`X`表示一个数组或者矩阵；`MARGIN`表示维度，1表示按行计算，2表示按列计算；`FUN`是需要施加在`X`的行或者列上的函数。

例如，我们希望得到矩阵`x`每一行的标准差。利用`for`循环，我们可以将函数写成

```r{.line-numbers}
myrowsd1 <- function(x){
  n <- nrow(x)
  p <- ncol(x)
  rowsd <- c()
  for(i in 1:n){
    rowsd <- c(rowsd, sqrt(sum((x[i,] - mean(x[i,]))^2)/(p - 1)))
  }
  return(rowsd)
}
```

有了前面的向量化意识后，我们可以试着将代码写成

```r{.line-numbers}
## 算法2
myrowsd2 <- function(x){
  n <- nrow(x)
  p <- ncol(x)
  rmean <- rowMeans(x)
  xbar <- matrix(rep(rmean, p), ncol = p)
  rowsd <- sqrt(rowSums((x - xbar)^2)/(p - 1))
  return(rowsd)
}
```

下面利用`apply`将`sd`函数作用到`x`的每一行上，立即得到我们想要的结果

```r{.line-numbers}
## 算法3
myrowsd3 <- function(x){
  return(apply(x, 1, sd))
}
```

从代码上看，`apply`函数最终的代码确实简洁漂亮。然而，是美与实力兼有还是花瓶一个，需要实践出真知。

为了比较这三种算法的差异，我们从正态分布中产生1000行10列的矩阵，重复100次后得到

```r
> x <- matrix(rnorm(10000), ncol = 10)
> system.time(# 算法1
+   replicate(100, myrowsd1(x))
+ )
用户 系统 流逝 
2.08 0.19 2.37 
> system.time(# 算法2
+   replicate(100, myrowsd2(x))
+ )
用户 系统 流逝 
0.03 0.04 0.08 
> system.time(# 算法3
+   replicate(100, myrowsd3(x))
+ )
用户 系统 流逝 
2.52 0.02 2.58 
```

结果令人大跌眼镜——`apply`函数耗时最多！使用内置函数`rowMeans`和`rowSums`的**算法2最高效。**

### - 并行apply函数

可否让`apply`函数的简洁与高效并存呢？考虑到`apply`使用`FUN`的过程是相互独立的，所以并行计算是个可行的思路。

在R中已经内置`parallel`程序包，可以使用类似`apply`但能多线程运算的`parRapply`（对矩阵的行进行处理）命令。它的基本调用形式为

```r
cl <- makeCluster(<size of pool>)# 建立并行线程数
parRapply(cl, x, FUN)
stopCluster(cl)# 关闭并行
```

我们利用并行`apply`命令重新计算上述问题得到

```r
> library(parallel)
> cl <- makeCluster(2)# 并行两次
> system.time(# 算法3
+   replicate(100, parRapply(cl, x, sd))
+ )
用户 系统 流逝 
0.39 0.16 2.33 
> stopCluster(cl)
```

可见并行两次的速度确实有提升，超过了`for`循环但远逊色于算法2的速度——我的surfaceGo只有两个物理核心，如果电脑允许你可以增加并行次数而降低运算时间。

**尽管如此，我想说的是`apply`函数并不是如我们以为的那样高效，我们在乎计算时间时便要警惕使用该命令！**

不过，这并不意味着`apply`函数没有存在的价值。在时间开销并不那么重要时，`apply`得到的代码十分优雅简洁，也容易理解。

### - apply族其他函数

`apply`族函数除了基本的`apply`命令外，常用的还有`lapply`和`sapply`等。

```r
lapply(X, FUN, ...)# X通常为向量、数据框或列表
```

`apply`命令并不适用于向量，而`lapply`则显得更加灵活——其中`X`可取向量、数据框、列表等。例如，我们给定一个向量，然后针对每个元素生成相应数目的均匀随机数，结果用列表储存。

```r
> x <- c(2,3,4)
> lapply(x, runif)
[[1]]
[1] 0.2485223 0.4080533

[[2]]
[1] 0.9620718 0.5607855 0.2601507

[[3]]
[1] 0.23525977 0.05190834 0.05177627 0.72261602
```

如果不喜欢列表储存结果，可以考虑`sapply`命令，它返回向量或者矩阵。如果`FUN`得到的结果是一维的，那么`sapply`就返回一个向量；如果`FUN`得到的结果是多维的，那么`sapply`就返回一个矩阵，其每一列对应一个输出结果。

例如，我们输入一个列表，求其中每一项的均值和方差。

```r
> x <- list(a = 1:3, b = rnorm(100))
> myfunction <- function(x){
+   return(c(mean(x), var(x)))
+ }
> sapply(x, myfunction)
     a          b
[1,] 2 -0.1352013
[2,] 1  1.1843918
```

更多的使用方式和注意细节限于篇幅，不再赘述。


## 3. 多线程运算

实际的数据分析中，我们常常需要重复某一过程数十次以上，例如Cross Validation寻找最优参数。一般每次计算相互独立，但计算的规律又极其相似。除了`for`循环外，我们能否有高效的方法实现该目的呢？

答案是肯定的：**R中的多线程运算就能办到。**

我们这里主要介绍`foreach`程序包中`foreach`命令的使用方式。假如我们需要计算100个向量的均值（计算互不影响），`for`循环意味着依次进行均值计算，这显然低效；`foreach`多线程则可以同时计算多个向量的均值（同时计算多少取决于你电脑的线程数），自然就显得高效。

因此，`foreach`可以看作`for`循环的加强版。它需要事先安装`foreach`和`doPrarllel`两个程序包。

```r{.line-numbers}
install.packages("foreach")
install.packages("doParallel")
```

### - 基本使用介绍

它的基本用法为

```r
library(foreach)
library(doParallel)
cl <- makeCluster(<size of pool>)
registerDoParallel(cl)# 注册并行的线程
foreach(..., .combine)  %dopar% FUN # 并行计算FUN
stopImplicitCluster()# 关闭并行
```

首先，我们需要确定多线程数并进行注册使用。如果你不知道自己电脑可用的线程数，可以使用`detectCores`命令获得；通过添加参数`logical = F`可以得到电脑的实际物理核心数，也即是真正可供调用的实体核心。

```r
> detectCores()
[1] 4
> detectCores(logical = F)
[1] 2
```

可见我的surfaceGo实际有两个物理核心，但可以虚拟两个核心出来。因此，我们可以设置的线程数不超过3个——电脑的运行也需要核心维持，故不能跑满。

其次，我们解释一下`foreach`的各个参数。`...`是函数`FUN`进行计算的变量，可为一个或多个。`.combine`则是设置输出的结果形式，如用`+`连接、用`rbind`拼接成矩阵等。

最后，我们使用`stopImplicitCluster()`关闭并行而不是`stopCluster()`命令。

### 案例分析

我们考虑两个向量对应的的分量和与分量积，输出结果为一个矩阵。那么，程序可以写成

```r{.line-numbers}
myfun <- function(a, b){
  return(c(a+b, a*b))
}
cl <- makeCluster(2)
registerDoParallel(cl)
res <- foreach(x = 1:5, y = 6:10, .combine = "rbind") %dopar% myfun(x, y)
stopImplicitCluster()
```

查看结果为

```r
> res
         [,1] [,2]
result.1    7    6
result.2    9   14
result.3   11   24
result.4   13   36
result.5   15   50
```

我在使用`foreach`包的过程中发现，**当问题规模不大、计算过程简单时，`foreach`调用计算机资源所消耗的时间反而让人得不偿失。** 我们用下面的案例来说明问题。

考虑前面的最小二乘案例，样本观测数固定为1000，重复估计1000次。

```r{.line-numbers}
> myfun <- function(){
+   x <- rnorm(1000)
+   y <- x*2 + 0.1*rnorm(1000)
+   ols3(x, y)
+ }
> system.time(# for循环
+   for(i in 1:1000){
+     myfun()
+   }
+ )
用户 系统 流逝 
0.36 0.00 0.36 
> cl <- makeCluster(2)
> registerDoParallel(cl)
> system.time(# foreach并行
+   foreach(i = 1:1000, .combine = "rbind") %dopar% myfun()
+ )
用户 系统 流逝 
0.83 0.33 1.45 
> stopImplicitCluster()
```

**这也和我们开头所呼应：如果算法先天不足、设计不合理，指望通过某个函数或程序包暴力提升运算速度，就很难得到满意的结果。**

### - 一点补充

现在我们编写函数，将`foreach`得到的100个最小二乘估计的平均值作为最终的解。你的函数若像下面一样编写

```r{.line-numbers}
meanb <- function(){
  cl <- makeCluster(2)
  registerDoParallel(cl)
  bhat <- foreach(i = 1:10, .combine = "rbind") %dopar% myfun()
  stopImplicitCluster()
  return(mean(bhat))
}
```

调用后，就现找不到`myfun`函数的错误提示

```r
> meanb()
 Error in myfun() : task 1 failed - "没有"myfun"这个函数"
```

**原因是`foreach`被嵌套进函数后无法使用函数外部的变量。** 因此，我们需要添加`.export`参数，将需要的外部变量传递到函数内部，也即是

```r{.line-numbers}
bhat <- foreach(i = 1:10, .export = c("ols3", "myfun"), .combine = "rbind")
```

这样，你在外部定义所需变量后再运行`meanb`函数，就能正常输出结果

```r
> meanb()
[1] 2.001765
}
```

## 4. 写在最后

本次我们学习了R提升速度的常见思路和注意事项。最后再次提醒大家，理解自己的问题和程序结构最重要，盲目通过某些命令来暴力提升计算速度，可能会适得其反。