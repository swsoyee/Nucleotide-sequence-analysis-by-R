解析 | 発現変動 | 2群間 | 対応なし | 複製あり | WAD(Kadota\_2008)
================
Kadota, Su
2018年5月4日

リファレンス：[解析 | 発現変動 | 2群間 | 対応なし | 複製あり |
WAD(Kadota\_2008)](http://www.iu.a.u-tokyo.ac.jp/~kadota/r_seq.html#analysis_deg_2_unpaired_ari_WAD)

[TCC](http://bioconductor.org/packages/release/bioc/html/TCC.html)パッケージを用いてiDEGES/edgeR正規化を行ったデータを入力として、
WAD法([Kadota et
al., 2008](http://www.ncbi.nlm.nih.gov/pubmed/18578891))を適用して発現変動遺伝子(Differentially
expressed Genes; DEGs)検出を行うやり方を示します。
WAD法は統計的手法ではない(ヒューリスティックな方法)ので、出力ファイル中にp-valueやq-valueは存在しません。
それゆえ、FDR閾値を満たす遺伝子数、という概念もありません。
「ファイル」−「ディレクトリの変更」で解析したいファイルを置いてあるディレクトリに移動し以下をコピペ。

## 1\. サンプルデータ13の10,000 genes × 6 samplesの[カウントデータ](http://www.iu.a.u-tokyo.ac.jp/~kadota/R_seq/data_hypodata_3vs3.txt)の場合：

-----

Biological replicatesを模倣したシミュレーションデータ(G1群3サンプル vs. G2群3サンプル)です。
gene\_1〜gene\_2000までがDEG (最初の1800個がG1群で高発現、残りの200個がG2群で高発現)
gene\_2001〜gene\_10000までがnon-DEGであることが既知です。

### 0\. 必要なパッケージをロード

``` r
libs <- c("knitr", "dplyr", "caret", "devtools", "plotly")
for (i in libs){
  if( !is.element(i, .packages(all.available = TRUE)) ) {
    install.packages(i)
  }
  library(i,character.only = TRUE)
}
library(TCC)
```

### 1\. 設置及びファイルの読み込み

``` r
in_f <- "http://www.iu.a.u-tokyo.ac.jp/~kadota/R_seq/data_hypodata_3vs3.txt"
# out_f1 <- "hoge1.txt"                  #出力ファイル名を指定してout_f1に格納
param_G1 <- 3                          #G1群のサンプル数を指定
param_G2 <- 3                          #G2群のサンプル数を指定
param_FDR <- 0.05                      #false discovery rate (FDR)閾値を指定

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

### 2\. 前処理（TCCオブジェクト作成）

``` r
#G1群を1、G2群を2としたベクトルdata.clを作成
data.cl <- c(rep(1, param_G1), rep(2, param_G2))
data.cl
```

    ## [1] 1 1 1 2 2 2

``` r
#TCCクラスオブジェクトtccを作成
tcc <- new("TCC", data, data.cl)
tcc
```

    ## Count:
    ##        G1_rep1 G1_rep2 G1_rep3 G2_rep1 G2_rep2 G2_rep3
    ## gene_1      36      56     144       2       1       0
    ## gene_2      84     152     124      52      37      28
    ## gene_3     592     840     800     151     257     200
    ## gene_4       0       8       4       1       1       3
    ## gene_5      32      32       0       1       1       0
    ## gene_6       4       0      24       4      10       0
    ## 
    ## Sample:
    ##         group norm.factors lib.sizes
    ## G1_rep1     1            1   1762346
    ## G1_rep2     1            1   1561258
    ## G1_rep3     1            1   1818047
    ## G2_rep1     2            1   1023545
    ## G2_rep2     2            1   1075566
    ## G2_rep3     2            1   1028008

``` r
#正規化を実行した結果をtccに格納
tcc <- calcNormFactors(tcc, 
                       norm.method="tmm", 
                       test.method="edger",
                       iteration=3, 
                       FDR=0.1, 
                       floorPDEG=0.05)
tcc
```

    ## Count:
    ##        G1_rep1 G1_rep2 G1_rep3 G2_rep1 G2_rep2 G2_rep3
    ## gene_1      36      56     144       2       1       0
    ## gene_2      84     152     124      52      37      28
    ## gene_3     592     840     800     151     257     200
    ## gene_4       0       8       4       1       1       3
    ## gene_5      32      32       0       1       1       0
    ## gene_6       4       0      24       4      10       0
    ## 
    ## Sample:
    ##         group norm.factors lib.sizes
    ## G1_rep1     1    0.7475138   1317378
    ## G1_rep2     1    0.8360228   1305247
    ## G1_rep3     1    0.7245591   1317282
    ## G2_rep1     2    1.2546831   1284225
    ## G2_rep2     2    1.1885056   1278316
    ## G2_rep3     2    1.2487157   1283690
    ## 
    ## DEGES:
    ##    Pipeline       : tmm - [ edger - tmm ] X 3
    ##    Execution time : 13.2 sec
    ##    Threshold type : FDR < 0.10
    ##    Potential PDEG : 0.13

### 3\. DEG検出

``` r
#DEG検出を実行した結果をtccに格納
tcc <- estimateDE(tcc, test.method="wad")
tcc
```

    ## Count:
    ##        G1_rep1 G1_rep2 G1_rep3 G2_rep1 G2_rep2 G2_rep3
    ## gene_1      36      56     144       2       1       0
    ## gene_2      84     152     124      52      37      28
    ## gene_3     592     840     800     151     257     200
    ## gene_4       0       8       4       1       1       3
    ## gene_5      32      32       0       1       1       0
    ## gene_6       4       0      24       4      10       0
    ## 
    ## Sample:
    ##         group norm.factors lib.sizes
    ## G1_rep1     1    0.7475138   1317378
    ## G1_rep2     1    0.8360228   1305247
    ## G1_rep3     1    0.7245591   1317282
    ## G2_rep1     2    1.2546831   1284225
    ## G2_rep2     2    1.1885056   1278316
    ## G2_rep3     2    1.2487157   1283690
    ## 
    ## DEGES:
    ##    Pipeline       : tmm - [ edger - tmm ] X 3
    ##    Execution time : 13.2 sec
    ##    Threshold type : FDR < 0.10
    ##    Potential PDEG : 0.13
    ## 
    ## Results:
    ##   gene_id  a.value    m.value p.value q.value rank estimatedDEG
    ## 1  gene_1 3.148234 -6.2619295      NA      NA  169            1
    ## 2  gene_2 6.096777 -1.5881518      NA      NA 1053            0
    ## 3  gene_3 8.601768 -1.8414834      NA      NA  329            1
    ## 4  gene_4 1.370464 -1.2335151      NA      NA 5335            0
    ## 5  gene_5 1.916709 -4.9665865      NA      NA 1701            0
    ## 6  gene_6 2.721477 -0.9585627      NA      NA 6489            0

``` r
#WAD統計量を抽出した結果をstatisticに格納
statistic <- tcc$stat$testStat
#WAD統計量でランキングした結果をrankingに格納
ranking <- tcc$stat$rank               
```

### 4\. ファイルに保存(テキストファイル)

``` r
#入力データの右側にDEG検出結果を結合したものをtmpに格納
tmp <- cbind(rownames(tcc$count), tcc$count, statistic, ranking)
kable(head(tmp))
```

|         |         | G1\_rep1 | G1\_rep2 | G1\_rep3 | G2\_rep1 | G2\_rep2 | G2\_rep3 | statistic            | ranking |
| ------- | ------- | :------- | :------- | :------- | :------- | :------- | :------- | :------------------- | :------ |
| gene\_1 | gene\_1 | 36       | 56       | 144      | 2        | 1        | 0        | \-1.1044909088107    | 169     |
| gene\_2 | gene\_2 | 84       | 152      | 124      | 52       | 37       | 28       | \-0.586546505782818  | 1053    |
| gene\_3 | gene\_3 | 592      | 840      | 800      | 151      | 257      | 200      | \-0.970899954123072  | 329     |
| gene\_4 | gene\_4 | 0        | 8        | 4        | 1        | 1        | 3        | \-0.0745189325186262 | 5335    |
| gene\_5 | gene\_5 | 32       | 32       | 0        | 1        | 1        | 0        | \-0.33637617513029   | 1701    |
| gene\_6 | gene\_6 | 4        | 0        | 24       | 4        | 10       | 0        | \-0.0476399290827648 | 6489    |

``` r
#tmpの中身を指定したファイル名で保存
# write.table(tmp, out_f, sep="\t", append=F, quote=F, row.names=F)
```

### 5\. AUC値を計算（Simulation Data Only）

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

    ## [1] 0.8419372

### 6\. WAD関数

iDEGES/edgeR正規化後のデータを明示的に取得して、
estimateDE関数ではなくTCCパッケージ中のWAD関数を用いてWAD法を実行するやり方です。
WAD法はlog2変換後のデータを入力とすることを前提としており、発現レベルに相当する数値が1未満のものを1に変換してからlogをとっています。

``` r
#正規化後のデータを取り出してnormalizedに格納
normalized <- getNormalizedData(tcc)
kable(head(normalized))
```

|         |  G1\_rep1 |   G1\_rep2 |   G1\_rep3 |   G2\_rep1 |   G2\_rep2 |   G2\_rep3 |
| ------- | --------: | ---------: | ---------: | ---------: | ---------: | ---------: |
| gene\_1 |  35.46198 |  55.675754 | 141.858196 |   2.020970 |   1.015156 |   0.000000 |
| gene\_2 |  82.74462 | 151.119903 | 122.155669 |  52.545218 |  37.560754 |  28.305369 |
| gene\_3 | 583.15254 | 835.136309 | 788.101089 | 152.583230 | 260.894964 | 202.181210 |
| gene\_4 |   0.00000 |   7.953679 |   3.940505 |   1.010485 |   1.015156 |   3.032718 |
| gene\_5 |  31.52176 |  31.814716 |   0.000000 |   1.010485 |   1.015156 |   0.000000 |
| gene\_6 |   3.94022 |   0.000000 |  23.643033 |   4.041940 |  10.151555 |   0.000000 |

``` r
#DEG検出を実行した結果をtccに格納
out <- WAD(normalized, data.cl, logged=F, floor=1)
kable(head(out))
```

|         |         wad | rank |
| ------- | ----------: | ---: |
| gene\_1 | \-1.1044909 |  169 |
| gene\_2 | \-0.5865465 | 1053 |
| gene\_3 | \-0.9709000 |  329 |
| gene\_4 | \-0.0745189 | 5335 |
| gene\_5 | \-0.3363762 | 1701 |
| gene\_6 | \-0.0476399 | 6489 |

``` r
#WAD統計量でランキングした結果をrankingに格納
ranking <- out$rank                    
```

-----

>   - WAD：[Kadota et al., Algorithms Mol.
>     Biol., 2008](http://www.ncbi.nlm.nih.gov/pubmed/18578891)
>   - [TCC](http://bioconductor.org/packages/release/bioc/html/TCC.html):
>     [Sun et al., BMC
>     Bioinformatics, 2013](http://www.ncbi.nlm.nih.gov/pubmed/23837715)
