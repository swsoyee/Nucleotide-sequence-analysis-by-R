解析 | 発現変動 | 2群間 | 対応なし | 複製あり | SAMseq(Li\_2013)
================
Kadota, Su
2018年5月2日

リファレンス：[解析 | 発現変動 | 2群間 | 対応なし | 複製あり |
SAMseq(Li\_2013)](www.iu.a.u-tokyo.ac.jp/~kadota/r_seq.html#analysis_deg_2_unpaired_ari_SAMseq)

[samr](http://cran.r-project.org/web/packages/samr/index.html)パッケージ中の[SAMseq
(Li and Tibshirani, 2013)](http://www-stat.stanford.edu/~tibs/SAM/)を用いて
発現変動遺伝子(Differentially expressed Genes; DEGs)検出を行うやり方を示します。
「ファイル」−「ディレクトリの変更」で解析したいファイルを置いてあるディレクトリに移動し以下をコピペ。

## 1\. サンプルデータ13の10,000 genes × 6 samplesの[カウントデータ](http://www.iu.a.u-tokyo.ac.jp/~kadota/R_seq/data_hypodata_3vs3.txt)の場合：

-----

Biological replicatesを模倣したシミュレーションデータ(G1群3サンプル vs. G2群3サンプル)です。
gene\_1〜gene\_2000までがDEG (最初の1800個がG1群で高発現、残りの200個がG2群で高発現)
gene\_2001〜gene\_10000までがnon-DEGであることが既知です。

### 0\. 必要なパッケージをロード

``` r
libs <- c("knitr", "dplyr", "caret", "devtools", "plotly", "samr")
for (i in libs){
  if( !is.element(i, .packages(all.available = TRUE)) ) {
    install.packages(i)
  }
  library(i,character.only = TRUE)
}
```

### 1\. 設置及びファイルの読み込み

``` r
in_f <- "http://www.iu.a.u-tokyo.ac.jp/~kadota/R_seq/data_hypodata_3vs3.txt"
# out_f1 <- "hoge1.txt"                  #出力ファイル名を指定してout_f1に格納
param_G1 <- 3                          #G1群のサンプル数を指定
param_G2 <- 3                          #G2群のサンプル数を指定

#入力ファイルの読み込み
data <- read.table(in_f, 
                   header=TRUE, 
                   row.names=1, 
                   sep="\t", 
                   quote="")           #in_fで指定したファイルの読み込み

kable(head(data, n=10))
```

|          | G1\_rep1 | G1\_rep2 | G1\_rep3 | G2\_rep1 | G2\_rep2 | G2\_rep3 |
| -------- | -------: | -------: | -------: | -------: | -------: | -------: |
| gene\_1  |       36 |       56 |      144 |        2 |        1 |        0 |
| gene\_2  |       84 |      152 |      124 |       52 |       37 |       28 |
| gene\_3  |      592 |      840 |      800 |      151 |      257 |      200 |
| gene\_4  |        0 |        8 |        4 |        1 |        1 |        3 |
| gene\_5  |       32 |       32 |        0 |        1 |        1 |        0 |
| gene\_6  |        4 |        0 |       24 |        4 |       10 |        0 |
| gene\_7  |      344 |      240 |      236 |       76 |       67 |       71 |
| gene\_8  |     1264 |      784 |     1060 |      212 |      183 |      179 |
| gene\_9  |       92 |       88 |       84 |       21 |       22 |       33 |
| gene\_10 |       64 |       48 |       96 |       24 |       13 |       12 |

### 2\. ラベル情報の作成

``` r
#G1群を1、G2群を2としたベクトルdata.clを作成
data.cl <- c(rep(1, param_G1), rep(2, param_G2))
data.cl
```

    ## [1] 1 1 1 2 2 2

### 3\. DEG検出

``` r
#DEG検出を実行した結果をoutに格納
out <- SAMseq(data, 
              data.cl, 
              resp.type="Two class unpaired")
```

    ## Estimating sequencing depths...
    ## Resampling to get new data matrices...
    ## perm= 1
    ## perm= 2
    ## perm= 3
    ## perm= 4
    ## perm= 5
    ## perm= 6
    ## perm= 7
    ## perm= 8
    ## perm= 9
    ## perm= 10
    ## perm= 11
    ## perm= 12
    ## perm= 13
    ## perm= 14
    ## perm= 15
    ## perm= 16
    ## perm= 17
    ## perm= 18
    ## perm= 19
    ## perm= 20
    ## perm= 21
    ## perm= 22
    ## perm= 23
    ## perm= 24
    ## perm= 25
    ## perm= 26
    ## perm= 27
    ## perm= 28
    ## perm= 29
    ## perm= 30
    ## perm= 31
    ## perm= 32
    ## perm= 33
    ## perm= 34
    ## perm= 35
    ## perm= 36
    ## perm= 37
    ## perm= 38
    ## perm= 39
    ## perm= 40
    ## perm= 41
    ## perm= 42
    ## perm= 43
    ## perm= 44
    ## perm= 45
    ## perm= 46
    ## perm= 47
    ## perm= 48
    ## perm= 49
    ## perm= 50
    ## perm= 51
    ## perm= 52
    ## perm= 53
    ## perm= 54
    ## perm= 55
    ## perm= 56
    ## perm= 57
    ## perm= 58
    ## perm= 59
    ## perm= 60
    ## perm= 61
    ## perm= 62
    ## perm= 63
    ## perm= 64
    ## perm= 65
    ## perm= 66
    ## perm= 67
    ## perm= 68
    ## perm= 69
    ## perm= 70
    ## perm= 71
    ## perm= 72
    ## perm= 73
    ## perm= 74
    ## perm= 75
    ## perm= 76
    ## perm= 77
    ## perm= 78
    ## perm= 79
    ## perm= 80
    ## perm= 81
    ## perm= 82
    ## perm= 83
    ## perm= 84
    ## perm= 85
    ## perm= 86
    ## perm= 87
    ## perm= 88
    ## perm= 89
    ## perm= 90
    ## perm= 91
    ## perm= 92
    ## perm= 93
    ## perm= 94
    ## perm= 95
    ## perm= 96
    ## perm= 97
    ## perm= 98
    ## perm= 99
    ## perm= 100
    ## Number of thresholds chosen (all possible thresholds) = 102
    ## Getting all the cutoffs for the thresholds...
    ## Getting number of false positives in the permutation...

``` r
#p値をp.valueに格納
p.value <- samr.pvalues.from.perms(out$samr.obj$tt,
                                   out$samr.obj$ttstar)  
#q値をq.valueに格納
q.value <- p.adjust(p.value, method="BH")  
#p.valueでランキングした結果をrankingに格納
ranking <- rank(p.value) 
#FDR閾値(q.value < 0.10)を満たす遺伝子数を表示
sum(q.value < 0.10)
```

    ## [1] 1030

### 4\. ファイルに保存(テキストファイル)

``` r
#入力データの右側にp.value、q.value、rankingを結合した結果をtmpに格納。
tmp <- cbind("gene_id"=rownames(data), data, p.value, q.value, ranking)

#tmpの中身を指定したファイル名で保存
# write.table(tmp, out_f, sep="\t", append=F, quote=F, row.names=F)
kable(head(tmp))
```

|         | gene\_id | G1\_rep1 | G1\_rep2 | G1\_rep3 | G2\_rep1 | G2\_rep2 | G2\_rep3 |   p.value |   q.value | ranking |
| ------- | :------- | -------: | -------: | -------: | -------: | -------: | -------: | --------: | --------: | ------: |
| gene\_1 | gene\_1  |       36 |       56 |      144 |        2 |        1 |        0 | 0.0099765 | 0.0968592 |   515.5 |
| gene\_2 | gene\_2  |       84 |      152 |      124 |       52 |       37 |       28 | 0.0099765 | 0.0968592 |   515.5 |
| gene\_3 | gene\_3  |      592 |      840 |      800 |      151 |      257 |      200 | 0.0099765 | 0.0968592 |   515.5 |
| gene\_4 | gene\_4  |        0 |        8 |        4 |        1 |        1 |        3 | 0.7189380 | 0.9446039 |  7544.5 |
| gene\_5 | gene\_5  |       32 |       32 |        0 |        1 |        1 |        0 | 0.2697250 | 0.7469538 |  3560.0 |
| gene\_6 | gene\_6  |        4 |        0 |       24 |        4 |       10 |        0 | 0.7897265 | 0.9536608 |  8197.0 |

### 5\. AUC値（と ConfusionMatrix）を計算（Simulation Data Only）

AUC値とConfusionMatrixを計算する。

``` r
# ROCパッケージはCRANにいないため、インストールされていない場合は下記のコードを実行してインストールしてください。
# source("https://bioconductor.org/biocLite.R")
# biocLite("ROC")
library(ROC)
param_DEG <- 1:2000                    #DEGの位置を指定
obj <- rep(0, nrow(data))           #初期値として全てが0の(non-DEGに相当)ベクトルobjを作成
obj[param_DEG] <- 1                    #DEGの位置に1を代入
AUC(rocdemo.sca(truth=obj, data=-ranking))#AUC計算
```

    ## [1] 0.8257229

``` r
param_FDR <- 0.1
deg_count <- sum(q.value < param_FDR)
prediction <- if_else(ranking <= deg_count, 1, 0)
confusionMatrix(prediction, obj)
```

    ## Confusion Matrix and Statistics
    ## 
    ##           Reference
    ## Prediction    0    1
    ##          0 7932 1038
    ##          1   68  962
    ##                                           
    ##                Accuracy : 0.8894          
    ##                  95% CI : (0.8831, 0.8955)
    ##     No Information Rate : 0.8             
    ##     P-Value [Acc > NIR] : < 2.2e-16       
    ##                                           
    ##                   Kappa : 0.5775          
    ##  Mcnemar's Test P-Value : < 2.2e-16       
    ##                                           
    ##             Sensitivity : 0.9915          
    ##             Specificity : 0.4810          
    ##          Pos Pred Value : 0.8843          
    ##          Neg Pred Value : 0.9340          
    ##              Prevalence : 0.8000          
    ##          Detection Rate : 0.7932          
    ##    Detection Prevalence : 0.8970          
    ##       Balanced Accuracy : 0.7363          
    ##                                           
    ##        'Positive' Class : 0               
    ## 

## 2\. サンプルデータ13の10,000 genes × 6 samplesの[カウントデータ](http://www.iu.a.u-tokyo.ac.jp/~kadota/R_seq/data_hypodata_3vs3.txt)の場合：

-----

「G1\_rep1, G1\_rep2, G1\_rep3, G2\_rep1, G2\_rep2,
G2\_rep3」の計6サンプル分からなります。
全10,000遺伝子中の最初の2,000個(gene\_1〜gene\_2000まで)が発現変動遺伝子(DEG)です。
全2,000
DEGsの内訳：最初の90%分(gene\_1〜gene\_1800)がG1群で4倍高発現、残りの10%分(gene\_1801〜gene\_2000)がG2群で4倍高発現
このDEG or non-DEGの位置情報と実際のランキング結果情報を用いてAUC値の計算を行う例です。

``` r
#必要なパッケージをロード
library(TCC)
library(samr)
library(ROC)

#シミュレーションデータの作成
set.seed(1000)
tcc <- simulateReadCounts(Ngene=10000, 
                          PDEG=0.2,
                          DEG.assign=c(0.9, 0.1), 
                          DEG.foldchange=c(4, 4),
                          replicates=c(3, 3))

#カウントデータ情報をdataに格納
data <- tcc$count
#G1群を1、G2群を2としたベクトルdata.clを作成
data.cl <- tcc$group$group             

#本番(DEG検出)
#DEG検出を実行した結果をoutに格納
out <- SAMseq(data, 
              data.cl, 
              resp.type="Two class unpaired")
```

    ## Estimating sequencing depths...
    ## Resampling to get new data matrices...
    ## perm= 1
    ## perm= 2
    ## perm= 3
    ## perm= 4
    ## perm= 5
    ## perm= 6
    ## perm= 7
    ## perm= 8
    ## perm= 9
    ## perm= 10
    ## perm= 11
    ## perm= 12
    ## perm= 13
    ## perm= 14
    ## perm= 15
    ## perm= 16
    ## perm= 17
    ## perm= 18
    ## perm= 19
    ## perm= 20
    ## perm= 21
    ## perm= 22
    ## perm= 23
    ## perm= 24
    ## perm= 25
    ## perm= 26
    ## perm= 27
    ## perm= 28
    ## perm= 29
    ## perm= 30
    ## perm= 31
    ## perm= 32
    ## perm= 33
    ## perm= 34
    ## perm= 35
    ## perm= 36
    ## perm= 37
    ## perm= 38
    ## perm= 39
    ## perm= 40
    ## perm= 41
    ## perm= 42
    ## perm= 43
    ## perm= 44
    ## perm= 45
    ## perm= 46
    ## perm= 47
    ## perm= 48
    ## perm= 49
    ## perm= 50
    ## perm= 51
    ## perm= 52
    ## perm= 53
    ## perm= 54
    ## perm= 55
    ## perm= 56
    ## perm= 57
    ## perm= 58
    ## perm= 59
    ## perm= 60
    ## perm= 61
    ## perm= 62
    ## perm= 63
    ## perm= 64
    ## perm= 65
    ## perm= 66
    ## perm= 67
    ## perm= 68
    ## perm= 69
    ## perm= 70
    ## perm= 71
    ## perm= 72
    ## perm= 73
    ## perm= 74
    ## perm= 75
    ## perm= 76
    ## perm= 77
    ## perm= 78
    ## perm= 79
    ## perm= 80
    ## perm= 81
    ## perm= 82
    ## perm= 83
    ## perm= 84
    ## perm= 85
    ## perm= 86
    ## perm= 87
    ## perm= 88
    ## perm= 89
    ## perm= 90
    ## perm= 91
    ## perm= 92
    ## perm= 93
    ## perm= 94
    ## perm= 95
    ## perm= 96
    ## perm= 97
    ## perm= 98
    ## perm= 99
    ## perm= 100
    ## Number of thresholds chosen (all possible thresholds) = 107
    ## Getting all the cutoffs for the thresholds...
    ## Getting number of false positives in the permutation...

``` r
#p値をp.valueに格納
p.value <- samr.pvalues.from.perms(out$samr.obj$tt,
                                   out$samr.obj$ttstar)
#q値をq.valueに格納
q.value <- p.adjust(p.value, method="BH")
#p.valueでランキングした結果をrankingに格納
ranking <- rank(p.value)
#FDR閾値(q.value < 0.10)を満たす遺伝子数を表示
sum(q.value < 0.10)
```

    ## [1] 1065

``` r
#ファイルに保存
#入力データの右側にp.value、q.value、rankingを結合した結果をtmpに格納。
tmp <- cbind(rownames(data), data, p.value, q.value, ranking)
#tmpの中身を指定したファイル名で保存
# write.table(tmp, out_f, sep="\t", append=F, quote=F, row.names=F)

#AUC値を計算(高いほどよい方法)
#DEGの位置を1、non-DEGの位置を0としたベクトルobjを作成
obj <- as.numeric(tcc$simulation$trueDEG != 0)
AUC(rocdemo.sca(truth=obj, data=-ranking))#AUC計算
```

    ## [1] 0.8622411

``` r
# CufusionMatrix
deg_count <- sum(q.value < 0.10)
prediction <- if_else(ranking <= deg_count, 1, 0)
confusionMatrix(prediction, obj)
```

    ## Confusion Matrix and Statistics
    ## 
    ##           Reference
    ## Prediction    0    1
    ##          0 7941  994
    ##          1   59 1006
    ##                                           
    ##                Accuracy : 0.8947          
    ##                  95% CI : (0.8885, 0.9007)
    ##     No Information Rate : 0.8             
    ##     P-Value [Acc > NIR] : < 2.2e-16       
    ##                                           
    ##                   Kappa : 0.601           
    ##  Mcnemar's Test P-Value : < 2.2e-16       
    ##                                           
    ##             Sensitivity : 0.9926          
    ##             Specificity : 0.5030          
    ##          Pos Pred Value : 0.8888          
    ##          Neg Pred Value : 0.9446          
    ##              Prevalence : 0.8000          
    ##          Detection Rate : 0.7941          
    ##    Detection Prevalence : 0.8935          
    ##       Balanced Accuracy : 0.7478          
    ##                                           
    ##        'Positive' Class : 0               
    ## 

SAMseq法は既にTCCパッケージに内装されていますので、下記のコードはTCCパッケージを用いています。

``` r
# #必要なパッケージをロード
# library(TCC)
# library(samr)
# library(ROC)
# 
# #シミュレーションデータの作成
# set.seed(1000)
# tcc <- simulateReadCounts(Ngene=10000, 
#                           PDEG=0.2,
#                           DEG.assign=c(0.9, 0.1), 
#                           DEG.foldchange=c(4, 4),
#                           replicates=c(3, 3))

#本番(DEG検出)
#DEG検出を実行した結果をtccに格納
tcc <- estimateDE(tcc, test.method="samseq")
tcc
```

    ## Count:
    ##        G1_rep1 G1_rep2 G1_rep3 G2_rep1 G2_rep2 G2_rep3
    ## gene_1     327     149      41      52      41      17
    ## gene_2      21      10      14       6       6       6
    ## gene_3      65      70     176       4       5      17
    ## gene_4     370     279     539      82      62     108
    ## gene_5      55      15      20       3      15      22
    ## gene_6     633     251     483      53     127      73
    ## 
    ## Sample:
    ##         group norm.factors lib.sizes
    ## G1_rep1     1            1   1565439
    ## G1_rep2     1            1   1713884
    ## G1_rep3     1            1   1629431
    ## G2_rep1     2            1   1044355
    ## G2_rep2     2            1   1082128
    ## G2_rep3     2            1   1059770
    ## 
    ## Results:
    ##   gene_id  a.value    m.value   p.value    q.value   rank estimatedDEG
    ## 1  gene_1 6.358574 -1.6298033 0.1346640 0.54431690 2435.5            0
    ## 2  gene_2 3.288355 -0.7158278 0.0549565 0.33798585 1611.5            0
    ## 3  gene_3 4.940584 -2.9584924 0.0097025 0.09101782  533.5            1
    ## 4  gene_4 7.549480 -1.6184978 0.0097025 0.09101782  533.5            1
    ## 5  gene_5 4.367127 -0.5829117 0.3012855 0.74594083 3981.0            0
    ## 6  gene_6 7.654766 -1.8379127 0.0097025 0.09101782  533.5            1

``` r
#p値などの計算結果をresultに格納
result <- getResult(tcc, sort=FALSE)
kable(head(result))
```

| gene\_id |  a.value |     m.value |   p.value |   q.value |   rank | estimatedDEG |
| :------- | -------: | ----------: | --------: | --------: | -----: | -----------: |
| gene\_1  | 6.358574 | \-1.6298033 | 0.1346640 | 0.5443169 | 2435.5 |            0 |
| gene\_2  | 3.288355 | \-0.7158278 | 0.0549565 | 0.3379859 | 1611.5 |            0 |
| gene\_3  | 4.940584 | \-2.9584924 | 0.0097025 | 0.0910178 |  533.5 |            1 |
| gene\_4  | 7.549480 | \-1.6184978 | 0.0097025 | 0.0910178 |  533.5 |            1 |
| gene\_5  | 4.367127 | \-0.5829117 | 0.3012855 | 0.7459408 | 3981.0 |            0 |
| gene\_6  | 7.654766 | \-1.8379127 | 0.0097025 | 0.0910178 |  533.5 |            1 |

``` r
#FDR閾値(q.value < 0.10)を満たす遺伝子数を表示
sum(q.value < 0.10)
```

    ## [1] 1065

``` r
#AUC値を計算(高いほどよい方法)
calcAUCValue(tcc)
```

    ## [1] 0.8583803

-----

>   - [SAMseq](http://www-stat.stanford.edu/~tibs/SAM/)：[Li and
>     Tibshirani, Stat Methods Med
>     Res., 2013](http://www.ncbi.nlm.nih.gov/pubmed/22127579)
>   - [TCC](http://bioconductor.org/packages/release/bioc/html/TCC.html):
>     [Sun et al., BMC
>     Bioinformatics, 2013](http://www.ncbi.nlm.nih.gov/pubmed/23837715)
>   - [samr](http://cran.r-project.org/web/packages/samr/index.html)
