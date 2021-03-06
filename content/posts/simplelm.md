---
title: "简单回归分析的R实现"
date: 2019-07-26
toc: true
categories:
  - 笔记
tags:
  - 多元统计
  - r语言
  - 回归分析
---

## 1 导言
回归分析是一个古老的话题。一百多年前，英国的统计学家高尔顿 (F. Galton，1822-1911) 和他的学生皮尔逊 (K. Pearson，1857-1936) 在研究父母和其子女身高的遗传关系问题中，统计了1078对夫妇的平均身高以及他们的一个成年儿子的身高数据。

他们将孩子的身高作为自变量 $x$，父母的平均身高作为因变量 $y$，然后将两者画在同一张直角坐标系上。结果，他们发现这些数据点“惊人的”位于一条直线的附近，并且经过计算得到了直线的拟合方程:

$$ y = 33.73 + 0.516x $$

> 这个结果看起来是违背直觉的。因为统计的结果表明，高个子父母的子女有低于父母身高的趋势；而矮个子的子女则有高于父母的趋势。高尔顿解释说，**自然界存在某种约束力将人的身高向某个“平均数”靠拢——或者说是“回归”——也即是统计学上回归的涵义。**

那么本文的主题便是了解线性回归模型并通过`R`来解决线性回归分析中的若干问题。


## 2 基础回顾
回归的概念来源于实际问题，那么现在我们所说的线性回归分析问题具体指的是什么呢？一般说来，如果我们研究的问题中的 $p$ 个自变量 $x_1$, $x_2$, ..., $x_p$ 和因变量 $y$ 的关系形式如下所示

$$ y_i = \beta_0 + \beta_1x_{i1} + \cdots + \beta_px_{ip} + \epsilon_i $$

那么我们就说这是一个线性回归问题，其中 $\epsilon_i$ 是随机误差项，$i$ 表示第 $i$ 个观测。在线性回归问题中我们的核心任务就是估计出未知参数 $\beta_0$, $\beta_1$, $\cdots$, $\beta_p$ 的值。

> 注意，线性回归问题的确定并不是通过自变量的形式，而是问题中待估计的未知参数最高次都为一次且关于未知参数呈线性关系。例如 $y = \beta_0 + \beta_1x_1^2 + \epsilon$；$y = \beta_0 + \beta_1x_1x_2 + \epsilon$ 都是线性回归问题。

通常在线性回归中估计未知参数方法是[最小二乘法（OLS）](/2019/08/forum-ols/)，而为了保证估计值能够很好的解释模型，我们又有如下前提条件：

- **正态性**：$\epsilon_i$ 服从正态分布；
- **独立性**：$\epsilon_i$ 之间是独立的；
- **线性性**：$x$ 和 $y$ 必须线性相关；
- **同方差性**：$\epsilon_i$ 的方差不变。

这些条件又被称为**高斯—马尔可夫条件，** 它们保证了在经典线性回归中最小二乘估计的优越性。


## 3 求解线性回归模型函数
### 3.1 极大似然法
最小二乘法和极大似然法都可以用来求解线性回归模型，我们在往期文章中讨论过最小二乘法，这里对似然法进行简单介绍。

假设我们得到下面一组观测数据：

$$ (x_1, y_1), (x_2, y_2), \cdots, (x_n, y_n) $$

那么根据高斯-马尔可夫假设，我们可以得到残差估计的**似然函数** 为

$$ \mathscr{L} = \prod_{i=1}^n \frac{1}{\sqrt{2\pi}\sigma}\exp\left[-\frac{(y_i - x_i^T\hat{\beta} - \hat{\beta}_0)^2}{2\sigma^2}\right] $$

> 这个式子的成立还需要假设残差分布的均值为0，标准差为 $\sigma$。这个假设是可行的。因为残差如果均值不为零，可以将其移到模型的截距项里。

如何通过上面的函数得到系数的估计值呢？**极大似然的思想便是，让这些估计值使得似然函数达到最大！** 这个想法很朴素：每个观测数据随机且互相独立，我们一次搜集便得到眼前的数据，那么**自然而然认为这些数据组合出现的概率是最大的。**

不过，数据已经搜集好便不能改动。**我们自然想到，系数的估计值便是让这些数据对应的概率可能性最大——也即是似然函数最大。**

现在假装大家已经理解了极大似然的原理，下面我们来求解它。直接最大化不太可行，我们通常对似然函数取对数得到对数似然函数

$$ \ln\mathscr{L} = -\frac{n}{2}\ln(2\pi\sigma^2) - \sum_{i = 1}^n\frac{(y_i - x_i^T\hat{\beta} - \hat{\beta}_0)^2}{2\sigma^2} $$

然后再分别对各个参数进行优化。限于篇幅，不再赘述。

### 3.2 R求解线性回归模型
我们可以利用现有软件进行模型求解。在R中求解线性回归问题的最基本函数就是`lm()`，其格式为：

```r{.line-numbers}
myfit <- lm(formula, data)
# formula 是要拟合的模型形式，用一个R公式表示
# data 就是模型的数据构成的数据框
```

下面我们解释一下`formula`具体的形式，首先看下表总结的`formula`中常用的符号

|符号 |说明 |
|:--- |:---|
|`~`  |分隔符号，左边为因变量，右边为自变量 |
|`+`  |分隔自变量 |
|`:`  |自变量的交互项，如 `xz` 可以表示成 `x:z` |
|`*`  |自变量的所有交互项，如 `x*z*w` 展开即为 `x+z+w+x:z+x:w+z:w+x:z:w` |
|`^`  |交互项可以达到某个次数，如`(x+z+w)^2`展开即为`x+z+w+x:z+x:w+z:w`|
|`.`  |除因变量外的所有自变量 |
|`-1`  |删除截距项 |
|`I()` |如`x+I((z+w)^2)`等价于`x+h`，`h`是`z+w`平方构成的新变量 |

如果自变量为 $x_1$, $x_2$ 和 $x_3$ 而预测变量为 $y$，我们假定的线性模型形式为：

$$ y = \beta_0 + \beta_1x_1 + \beta_2x_2 + \beta_3x_3 + \epsilon $$

那么`formula`可以写成：

```r{.line-numbers}
y ~ x1 + x2 + x3
# 或者为 y ~ .
```

其他形式的模型`formula`的表达式还请读者自行琢磨。

当模型拟合成功后，我们使用`summary()`函数来得到拟合的具体结果。而其他常用的获取线性回归模型拟合结果的函数如下表所示。

|函数 |说明 |
|:--- |:--- |
|`summary()` |拟合情况概述，包括系数显著性、模型显著性等 |
|`coefficients()` |拟合系数 |
|`fitted()` |因变量拟合值 |
|`residuals()` |残差值 |
|`plot()` |画出模型诊断图 |

## 4 实例分析
下面我们将用实例具体介绍`lm()`函数的使用方法。

### 4.1 简单线性回归
本例中我们使用基础安装中的数据集`women`数据，它记录了15个年龄在30~39岁间女性的身高和体重信息，我们现在来探究体重关于身高的关系。

```r{.line-numbers}
myfit <- lm(weight~height, data = women)
summary(myfit) # 展示拟合详细结果
```

程序的输出结果如下所示

![](/lm/simpleols.png)

这里主要给读者解释这么几项指标的含义：

- `Residuals`  体重预测值和真实值之差的统计信息，从左到右分别为最小值、下四分位数、中位数、上四分位数和最大值。

- `Coefficients` 第一列`Estimate`中`Intercept`对应的数值为截距项，height对应的即为身高变量前的估计系数。

- `Multiple R-squared`  介于0-1之间，越接近1说明线性关系越强。

- `p-value`  模型的F检验统计量的`p`值，值越小说明模型越可靠。

因此本例中体重和身高的回归方程为：

$$ \hat{Weight} = -87.51667 + 3.45\times Height $$

根据`R`方 (`Multiple R-squared`) 和`p`值 (`p-value`) 可知模型是可靠的。此外，我们可以作图观察最终的拟合结果。

![](/lm/simplefig.jpg)

### 4.2 具有交互项的线性回归
继续考虑上例，如果模型中存在一个交互项比如一个平方项，那么即有：

```r{.line-numbers}
myfit <- lm(weight~height + I(height^2), data = women)
summary(myfit) # 展示拟合详细结果
```

程序的输出结果如下所示。

![](/lm/jiaohu.png)

可以看到通过比较`R`方、`p`值，添加了平方项的线性模型效果更好。我们同样可以做出相应的图像。
![](/lm/jiaohufig.jpg)


## 5 写在最后
本文主要介绍了`R`中线性回归分析的简单操作方法。不过，这里仅仅涉及线性回归分析的冰山一角，关于线性回归问题中的回归诊断和异常点的判断等内容，限于篇幅这里就不做介绍了。有兴趣的读者可以学习《R in action》第8章中关于回归的讲解。
