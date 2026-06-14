# Lambda 选择说明

本项目使用的 nuclear-norm estimator 形式为：

```text
min_M 0.5 * mean((y - X vec(M))^2) + lambda * nuclear_norm(M)
```

代码中的默认 lambda 为：

```text
lambda = lambda_factor * sigma * sqrt((p + q) / n)
```

默认 `lambda_factor=1` 在 `p=30,q=25,n=500,cov-scale=0.2` 下对应：

```text
lambda = 0.16583124
```

这个值在当前设计强度下偏大。直接后果是：

- singular values 被过度 shrinkage；
- final rank 下降到 1；
- true Frobenius error 停在 7 附近。

## 单 seed 现象

在 `seed=20260527` 下：

| lambda_factor | full F-error | rank-2 F-error | rank |
|---:|---:|---:|---:|
| 1.000 | 7.027 | 7.027 | 1 |
| 0.050 | 4.343 | 4.110 | 9 |
| 0.028 | 4.791 | 3.891 | 12 |

调小 lambda 后，full estimator 的 bias 明显下降，但会引入更多非零 singular values。对 rank 已知的 simulation，后接 rank-2 truncation 可以显著改善误差。

## 多 seed 推荐

30 个 seed 的经验推荐是：

| 用法 | lambda_factor | lambda | median F-error | q75 F-error |
|---|---:|---:|---:|---:|
| full nuclear estimator | 0.04 | 0.006633 | 4.456 | 4.723 |
| rank-2 truncation | 0.03 | 0.004975 | 4.025 | 4.334 |

解释上，`0.04` 是 full nuclear 的稳健折中；`0.03` 更适合后续做 rank-2 truncation，因为它保留了更多信号方向，之后由 truncation 去掉较小的噪声方向。

## 不应过度解释

这些值不是 universal constants。当前结论依赖：

```text
p = 30
q = 25
n = 500
cov-scale = 0.2
sigma = 0.5
rank = 2
```

如果 `cov-scale`、样本量或维度改变，推荐 `lambda_factor` 也可能改变。
