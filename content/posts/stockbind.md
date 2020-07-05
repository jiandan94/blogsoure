---
title: "股票数据的合并"
date: 2020-03-27
toc: true
categories:
  - R语言
tags:
  - R语言
  - 算法
---


## 1 背景介绍
最近，有朋友跟我说起股票数据处理的问题。由于课程要求，学生需在软件中下载**股票指数数据** 和对应的**成分股收盘价数据**，然后对指数和收盘价做一个分析，比如回归分析。

金融软件导出的数据是指数+所有的成分股的CSV文件，每个CSV文件格式相同，包含时间、开盘价、收盘价、成交量等。那么，问题来了——怎么把这些文件整合到一起去呢？

准确来说，这个问题要考虑以下这些点：

1. 指数和所有成分股要按时间一一对应，每一行代表一个时间点的指数和成分股收盘价数据；
2. 成分股遇到停盘数据便会有缺失，要能够根据最近的收盘价进行填补。

以前也有相应的R处理代码供使用。不过在和同学交流时发现，这个代码运行经常报错。因此，我自己就编了一个。下面说说我编写的思路。

## 2 算法设计
我设计这个算法的思路很朴素。首先，指数每天都有不会缺失。因此，我首先把指数都读进来提取出**指数的时间点** 和**收盘信息。** 然后，根据成分股的个数创建合适的矩阵先保存好指数数据。如下图所示。

![](/stock-organization/index.png)

接着，对第一个成分股进行处理——把**时间** 和**收盘价** 提取出来。先不用管它哪里缺失了需不需要填充数据。

> 我们要明确一点，导出的时间段内，指数肯定是最完整的。因此，我们将成分股按照时间跟它匹配肯定可以实现。

那么，我们就根据**时间点** 把第一个成分股匹配到指数所在的矩阵中，也就是下面这样的图像。

![](/stock-organization/stock1.png)

从上图看出，第一支股票在时间点3和时间点4有缺失数据。那么按照填充要求，最近的时间点2有数据，因此可以用来填充，就得到了下面的结果。

![](/stock-organization/addstock1.png)

以此类推，我们把有的股票都按照时间匹配到这个矩阵中，并做好缺失数据填充，结果就是下面这样。

![](/stock-organization/addstockp.png)

那么最后就简单了——按照**木桶原理** 把前面没有信息的数据全部砍掉即可，如下图所示。

![](/stock-organization/allinone.png)


## 3 R 程序代码
这个程序我当时为了赶作业图省事用了很多`for`循环，因此数据量大的时候，运行有点慢。

程序分为**日线数据处理代码** 和**分钟线数据处理代码。** 显然，可以把这两个合并到一起，不过课已上完我就懒得弄了。

另外，这是根据**西南金点子软件** 导出数据结构编写的，其他金融软件导出的数据对应改下应该就没问题了。

**- 代码使用方式**

```r{.line-numbers}
【1】我桌面是默认的工作目录
【2】上证指数日线csv数据放在桌面上
【3】桌面上的data文件夹用来放置50只成分股csv数据（一定要保存成csv格式）

那么参数就可以写成
wdir<-"C:/Users/***/Desktop/" ######## (要换成自己的路径)
filedir <- "data"
index.name <- "sz50.csv"

data_org(wdir,filedir,index.name)## 日线数据函数
minutedata_org(wdir,filedir,index.name)## 分钟线数据函数

其中出现下面警告直接无视，这只是将字符转换成数字时系统的温馨提示。
Warning message:
In data_org("C:/Users/tom/Desktop/", "data", "sz50.csv") :
  NAs introduced by coercion
  
程序运行后就将分散的指数和成分股数据按照时间整合到一起了，其中缺失数据用前一天补充。
```

**- 日线处理程序**

```r{.line-numbers}
data_org <- function(wdir,filedir,index.name){
  
  #设定文件所在路径与读入文件
  setwd(wdir)
  f.names <- list.files(filedir)
  stock_names <- substr(f.names,4,9)
  f.dir <- paste("./",filedir,"/",f.names,sep="")
  n.names <- length(f.names)
  
  # 首先阅读股票指数数据
  index_data <- read.csv(index.name,skip = 2,header = F,stringsAsFactors = F)
  index_data <- index_data[1:(nrow(index_data) - 1),c(1,5)]
  colnames(index_data) <- c("date","index")
  index_data$date <- as.Date(index_data$date)# 转换成时间格式
  index_data$index <- as.numeric(index_data$index)# 将指数数据转换成数值型
  index_data <- na.omit(index_data)# 去除NA的行
  
  # 处理成分股数据
  n_zero <- rep(0,n.names)
  stocks <- list()
  file_data <- index_data
  for(i in 1:n.names){
    # 读取数据
    stocks[[i]] <- read.csv(file = f.dir[[i]],skip = 2,header = F,stringsAsFactors = F)
    stock_data <- stocks[[i]][1:(nrow(stocks[[i]]) - 1),c(1,5)]
    colnames(stock_data) <- c("date",paste0("Stock",stock_names[i]))
    stock_data$date <- as.Date(stock_data$date)# 转换成时间格式
    stock_data[,2] <- as.numeric(stock_data[,2])# 将股票数据转换成数值型
    stock_data <- na.omit(stock_data)# 去除NA的行
    
    # 根据指数的日期来处理成分股数据
    temp1 <- rep(0,length(index_data$date))
    for (t in 1:length(stock_data$date)){
      date_index <- grep(stock_data$date[t],index_data$date)
      temp1[date_index] <- stock_data[,2][t]
    }
    
    for (k in 2:length(index_data$date)) {
      if(temp1[k] == 0){
        temp1[k] <- temp1[k - 1]
      }
    }
    
    n_zero[i] <- which(temp1 > 0)[1]
    temp2 <- data.frame(temp1)
    colnames(temp2) <- paste0("Stock",stock_names[i])

    file_data <- cbind(file_data,temp2)
  }
  
  data_trunc <- max(n_zero)
  file_data <- file_data[-(1:data_trunc),]
  
  write.csv(file_data,paste0("org_",index.name))
  
}
```

**- 分钟线处理程序**
```r{.line-numbers}
minutedata_org <- function(wdir,filedir,index.name){
  
  #设定文件所在路径与读入文件
  setwd(wdir)
  f.names <- list.files(filedir)
  stock_names <- substr(f.names,4,9)
  f.dir <- paste("./",filedir,"/",f.names,sep="")
  n.names <- length(f.names)
  
  # 首先阅读股票指数数据
  index_data <- read.csv(index.name,skip = 2,header = F,stringsAsFactors = F)
  index_data <- index_data[1:(nrow(index_data) - 1),c(1,2,6)]
  colnames(index_data) <- c("date","minute","index")
  index_data$dateminute <- paste0(index_data$date,index_data$minute)
  index_data$index <- as.numeric(index_data$index)# 将指数数据转换成数值型
  index_data <- na.omit(index_data)# 去除NA的行
  
  # 处理成分股数据
  n_zero <- rep(0,n.names)
  stocks <- list()
  file_data <- index_data
  for(i in 1:n.names){
    # 读取数据
    stocks[[i]] <- read.csv(file = f.dir[[i]],skip = 2,header = F,stringsAsFactors = F)
    stock_data <- stocks[[i]][1:(nrow(stocks[[i]]) - 1),c(1,2,6)]
    colnames(stock_data) <- c("date","minute",paste0("Stock",stock_names[i]))
    stock_data$dateminute <- paste0(stock_data$date,stock_data$minute)
    stock_data[,3] <- as.numeric(stock_data[,3])# 将股票数据转换成数值型
    stock_data <- na.omit(stock_data)# 去除NA的行
    
    # 根据指数的日期来处理成分股数据
    temp1 <- rep(0,length(index_data$dateminute))
    for (t in 1:length(stock_data$dateminute)){
      date_index <- grep(stock_data$dateminute[t],index_data$dateminute)
      temp1[date_index] <- stock_data[,3][t]
    }
    
    for (k in 2:length(index_data$dateminute)) {
      if(temp1[k] == 0){
        temp1[k] <- temp1[k - 1]
      }
    }
    
    n_zero[i] <- which(temp1 > 0)[1]
    temp2 <- data.frame(temp1)
    colnames(temp2) <- paste0("Stock",stock_names[i])
    
    file_data <- cbind(file_data,temp2)
  }
  
  data_trunc <- max(n_zero)
  file_data <- file_data[-(1:data_trunc),]
  
  write.csv(file_data,paste0("org_",index.name))
  
}
```