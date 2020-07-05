---
title: "漫谈最小二乘法"
date: 2019-08-02
toc: true
categories:
  - 学术杂谈
tags:
  - 数学史
---

> 如果随便逮到一个统计专业的学生问他“统计方法谁家强”，相信大部分人会异口同声得说“最小二乘法”。

的确，最小二乘法是一种非常重要的统计方法，它的重要性不仅仅体现在对问题求解的自然、简单、有效层面，其背后所蕴含的“最小二乘思想”更在不同领域、不同问题中应用广泛。

虽然以现在的眼光来看最小二乘法出发点朴素而又自然，但是它的产生也是历经波折。本文我们一起来了解一下最小二乘法的“心路历程”吧！

## 1 天文和测地学
早期，在天文和测地学中经常会遇到这么一种数据分析情况：我们有若干个可以测量的量 $x_1$, $x_2$,..., $x_p$ 和 $y$，他们之间呈现一种线性关系

$$ y = \beta_0 + \beta_1x_1 + \cdots + \beta_px_p $$

这里 $\beta_1$,..., $\beta_p$ 都是未知参数，需要我们运用一定的方法估计出来从而应用到实际问题中去。

![](/forum-ols/tianwen.jpg)

现在稍有统计背景的都知道这用最小二乘可以轻易求解，不过当时还没有最小二乘的概念呢，因此求解就是一个令人头疼的问题。

> 可能有人说根据线性方程理论这也不是一个很困难的问题啊！只要测量得到 $p$ 组数据，那么应用线性方程的知识不就将未知参数求解出来了么？

可是在实际天文和测地学中同一问题研究人员会测量多组数据以降低测量过程中产生的误差，而相应的数据量基本都是大于未知参数的个数的 (其实就是超定线性方程组)，这就使得问题的求解比较棘手。

我们希望的是，一方面对于方程的个数大于未知参数个数的问题，求解方法应该尽可能的简单、有效；另一方面，由于实际问题中测量误差不可避免，求解方法在理论上应该对误差控制有一个保证。所以总结起来就是四个字：**又快又稳！**

## 2 早期工作
在最小二乘还未诞生的时代，各路英雄豪杰面对该问题也是使劲浑身解数。不过他们的方法核心思路是从大量的方程中挑选或者组合恰当数量的方程来进行参数求解。**因此，各家各派使得都是挑选和组合的功夫。**

例如在1750年，天文学家梅耶发表了一种测定航海船只经度的方法，其中要从27个方程中求解出3个未知参数。**他主要采用分组法将27个方程分成3组分别相加得到3个方程，最终求解出未知参数。** 这个方法曾经一度流行，被冠以梅耶的名字。

![](/forum-ols/oula.jpg)

> 此外，拉普拉斯和欧拉也在天文学中研究过类似的问题，但是**解法极其繁杂而且杂乱无章。** 他们二人这样的大数学家一生不知道解决过多少数学中的“疑难杂症”，对于这么一个并不是很困难的问题竟然束手无策。这确实让人感到不可思议。

## 3 勒让德和最小二乘
勒让德是法国大数学家，在数学的很多领域包括椭圆几分、数论和几何等方面都有着重大的贡献。

$$ (\hat{\beta}, \hat{\beta}_0) = \min \sum _{i=1}^n \left( y_i - x_i^T\hat{\beta} - \hat{\beta}_0 \right)^2 $$

最小二乘法最早于1805年勒让德公开发表的文章《计算彗星轨道的新方法》中问世。在这本著作的附录中，勒让德描述了最小二乘法的思想、具体做法和优缺点。

> 最小二乘法一经提出，由于其思想自然合理、操作简单有效，很快就得到欧洲一些国家的天文和测地工作者的广泛使用。据不完全统计，自1805年到1864年，有关这一方法的研究论文约有250篇。

尽管勒让德的工作没有涉及到最小二乘的误差分析理论，**但是他也注意到了各个方程因为误差不独立而不能直接运用最小二乘法，** 这的确难能可贵。最小二乘法的“快”勒让德已经说明了，关于它的“稳”则是高斯的工作了。

## 4 高斯的工作
1809年，高斯在《绕日天体的运动理论》的一节中讨论了“数据结合”的问题，实际就是误差分布的确立问题。假设真值为 $\theta$，有 $n$ 个独立的测量值 $x_1$,..., $x_n$，高斯将后者的概率定为

$$ \mathscr{L}(\theta) = f(x_1 - \theta)\cdots f(x_n - \theta) $$

其中 $f$ 就是待定的误差密度函数。**在确立密度函数形式过程中，高斯有两个创新点。**

- 一是他没有采取已有的贝叶斯推理方法，而是直接将 $\mathscr{L}(\theta)$ 的最大值——**极大似然的思想** ——定为 $\theta$ 的估计值。

- 二是他先承认了观测值 $x_1$,..., $x_n$ 的算数平均值为 $θ$ 的估计值，然后再去找误差的密度函数来迎合这一点——在这样的 $f$ 下，$\theta$ 的估计值就是算数平均值。最后他得出只有在

$$ f(x) = \frac{1}{\sqrt{2\pi}h}e^{-\frac{1}{2h^2}} $$

时才成立。这就是均值为 $\theta$ 标准差为 $h$ 的正态分布.

> 使用正态分布就可以对最小二乘给出一种解释，也就是可以对其误差做出理论上的分析，保证了这种方法的优越性。后世将最小二乘的发明权归功于它，也正是因为这一项工作。

尽管高斯讨论最小二乘法的文章发表较晚，但是他声称自己很早之前就运用勒让德的最小二乘法来解决问题。这也导致两人后来最小二乘的首创权争论。

![](/forum-ols/guass.jpg)

不过在高斯的证明中有点**循环论证** 的感觉，先承认算数平均值估计的优越性，再得到误差正态密度函数形式，然后再说明算术平均值作为估计的合理性。这一缺陷在**拉普拉斯运用其发现的中心极限定理得以解决。** 他指出现实中的误差可以看成很多量的叠加，那么根据他的中心极限定理，误差的分布就是正态分布。

## 5 写在最后
从最小二乘法的发展历史来看，**一项科学理论的发展并无坦途，** 尽管这项理论看起来朴素而又简单。
