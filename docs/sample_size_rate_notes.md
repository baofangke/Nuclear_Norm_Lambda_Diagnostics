# Sample Size Rate 诊断

这个实验固定

```text
p = 30
q = 25
rank = 2
sigma = 0.5
signal = 6.0,4.4
design = kronecker
cov-scale = 0.2
seed = 20260527
```

只改变样本量

```text
n = 300, 400, 500, 700, 1000, 1500, 2000
```

为了减少不同样本量之间的随机波动，脚本先生成 `max(n)=2000` 的同一个数据集，然后对每个样本量使用前 `n` 个样本。lambda 因子不随 `n` 重新调参，而是沿用 `n=500` 诊断得到的经验选择：

```text
full nuclear estimator: lambda_factor = 0.04
rank-2 truncation:     lambda_factor = 0.03
```

每个样本量下仍使用

```text
lambda = lambda_factor * sigma * sqrt((p + q) / n)
```

## 结果

| n | full iter | full F-error | full adjusted ratio | rank iter | rank-2 F-error | rank adjusted ratio |
|---:|---:|---:|---:|---:|---:|---:|
| 300 | 222 | 4.959 | 3.28 | 261 | 4.654 | 3.07 |
| 400 | 386 | 4.844 | 3.70 | 613 | 4.528 | 3.45 |
| 500 | 206 | 4.403 | 3.75 | 620 | 4.178 | 3.56 |
| 700 | 360 | 4.079 | 4.12 | 709 | 3.629 | 3.66 |
| 1000 | 362 | 3.476 | 4.19 | 816 | 3.056 | 3.69 |
| 1500 | 221 | 2.747 | 4.06 | 434 | 2.483 | 3.67 |
| 2000 | 652 | 2.410 | 4.11 | 887 | 1.935 | 3.30 |

其中 adjusted ratio 使用

```text
adjusted minimax rate = sigma * sqrt(rank * (p + q) / n) / cov_scale
```

对 `log(||Mhat-M0||_F)` 关于 `log(n)` 做线性拟合：

```text
full nuclear slope       = -0.396
rank-2 truncation slope  = -0.464
1/sqrt(n) reference      = -0.5
```

![nuclear sample-size rate](../figures/fig_nuclear_sample_size_rate.png)

## 解释

在这组单 seed 结果中，两个 nuclear 估计量的 Frobenius error 都随样本量下降。rank-2 truncation 的误差整体更低，log-log 斜率也更接近 `-1/2`，说明在 rank 已知的 simulation 场景中，先用较小 lambda 保留信号方向，再做 rank-2 truncation 可以更接近 `1/sqrt(n)` 的经验下降关系。

full nuclear estimator 的下降斜率较浅，主要反映 nuclear penalty shrinkage 与额外非零 singular values 的折中。这个实验没有对每个 `n` 重新选择 lambda factor，所以它检验的是固定经验规则的 rate 行为，而不是每个样本量下的 oracle-tuned 最优表现。

更稳健的 rate 判断仍需要对每个样本量运行多个 seed，再汇总平均值和误差条。
