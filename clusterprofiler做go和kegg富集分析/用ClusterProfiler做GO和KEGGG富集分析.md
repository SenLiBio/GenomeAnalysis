# 用ClusterProfiler做GO和KEGGG富集分析

## 0. 一些参考教程

1.   李玉龙的一个总结：https://www.yuque.com/docs/share/e41471eb-5b6f-4c6a-90e2-62382c6818f0?#（密码：ykqg）
2.   徐州更的github：https://github.com/xuzhougeng/Learn-Bioinformatics/blob/master/8.Enrichment-Analysis.md
3.   Y叔公众号：http://guangchuangyu.github.io/2016/01/go-analysis-using-clusterprofiler/
4.   Y叔个人blog：https://guangchuangyu.github.io/cn/categories/
5.   一篇非常棒的blog：https://cloud.tencent.com/developer/article/2112467?from=article.detail.1457254
6.   非模式生物做富集：https://mp.weixin.qq.com/s/Mr3YLoc_-Y1WeLKJku1TzQ

## 1. 安装

```R
# 安装biomanager
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(version = "3.15")

# 安装clusterProfiler
BiocManager::install("clusterProfiler")

# 如果遇到有包报错没有装好，就手动装下，比如报错：不存在叫“DBI”这个名字的程序包
BiocManager::install(“DBI”)

# 安装完成
library("clusterProfiler")
```

## 2.模式生物做富集（或者说clusterprofiler已有的物种）

[note]: 有db数据库的物种列表可以从这个网站获得（关键词org）：https://bioconductor.org/packages/3.5/data/annotation/
[note2]: 这里拟南芥为例。

```bash
# 安装OrgDb文件
BiocManager::install("org.At.tair.db")

# 读入genelist，这里就是每行一个基因，不需要表头
gene_list<-read.table(file="Ath.genelist.txt",header=F)

# 富集分析
library(clusterProfiler)
ath_cc<-enrichGO(gene=gene_list$V1,OrgDb = "org.At.tair.db", keyType = "TAIR",ont="CC",pAdjustMethod = "BH", pvalueCutoff = 0.01,qvalueCutoff = 0.05)
ath_bp<-enrichGO(gene=gene_list$V1,OrgDb = "org.At.tair.db", keyType = "TAIR",ont="BP",pAdjustMethod = "BH", pvalueCutoff = 0.01,qvalueCutoff = 0.05)
ath_mf<-enrichGO(gene=gene_list$V1,OrgDb = "org.At.tair.db", keyType = "TAIR",ont="MF",pAdjustMethod = "BH", pvalueCutoff = 0.01,qvalueCutoff = 0.05)

# kegg分析
bruguiera_kegg<-enrichKEGG(gene=gene_list$V1,organism = 'ath', keyType = "kegg",pAdjustMethod = "BH", pvalueCutoff = 0.01,qvalueCutoff = 0.05)

# 画图，这里是一个dotplot的例子，具体的可以参考上面的那些教程
dotplot(bruguiera_kegg)
```

## 3.非模式生物做富集

[note]: 主要可以参考y叔的这篇笔记：https://mp.weixin.qq.com/s/Mr3YLoc_-Y1WeLKJku1TzQ

### 3.1 文件格式

#### 3.1.1 `go.db`文件

[note]: 该文件用文件夹下的做好的文件即可

| GO         | Description                                              | level              |
| ---------- | -------------------------------------------------------- | ------------------ |
| GO:0000001 | mitochondrion inheritance                                | biological_process |
| GO:0000007 | low-affinity zinc ion transmembrane transporter activity | molecular_function |
| ...        | ...                                                      | ...                |

#### 3.1.2 `GOannotation.tsv`

| Gene       | GO         | level              |
| ---------- | ---------- | ------------------ |
| Pg_S3686.2 | GO:0000165 | biological_process |
| Pg_S3686.2 | GO:0003674 | molecular_function |
| ...        | ...        | ...                |

#### 3.1.3 `KOannotation.tsv`

| Gene       | KO     | pathway  | decription         |
| ---------- | ------ | -------- | ------------------ |
| Pg_S6540.1 | K05907 | map00920 | Sulfur metabolism  |
| Pg_S6540.1 | K05907 | map01100 | Metabolic pathways |
| …          | …      | …        | …                  |

### 3.2 GO和kegg富集

==记得看y叔的这篇：https://mp.weixin.qq.com/s/Mr3YLoc_-Y1WeLKJku1TzQ==

```R
library(clusterProfiler)
KOannotation <- read.delim("KOannotation.tsv", stringsAsFactors=FALSE)
GOannotation <- read.delim("GOannotation.tsv", stringsAsFactors=FALSE)
GOinfo <- read.delim("go.tb", stringsAsFactors=FALSE)
# 前面获取gene list的过程略
gene_list<- # 你的gene list
# GO富集
## 拆分成BP，MF，CC三个数据框
GOannotation = split(GOannotation, with(GOannotation, level))
## 以MF为例
enricher（gene_list,
          TERM2GENE=GOannotation[['molecular_function']][c(2,1)],
          TERM2NAME=GOinfo[1:2]）
# KEGG富集
enricher（gene_list,
          TERM2GENE=KOannotation[c(3,1)],
          TERM2NAME=KOannotation[c(3,4)])
```

