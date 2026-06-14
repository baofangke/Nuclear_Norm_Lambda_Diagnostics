# Objective 收敛与真实 Frobenius Error

nuclear-norm estimator 的优化目标是 penalized objective：

```text
0.5 * mean((y - X vec(M))^2) + lambda * nuclear_norm(M)
```

simulation 中关心的诊断量通常是：

```text
Frobenius error to truth = norm(M_k - M0, "F")
```

这两个量不是同一个目标。因此，objective 持续下降或 objective-change 变小，并不保证真实 Frobenius error 也在同一时刻最小。

## 当前实验中的现象

在 `seed=20260527` 下，推荐 lambda 的迭代轨迹显示：

| lambda_factor | estimator | min F-error | min iter | final F-error | stop iter |
|---:|---|---:|---:|---:|---:|
| 0.04 | full nuclear | 4.342 | 267 | 4.494 | 427 |
| 0.04 | rank-2 truncation | 3.922 | 259 | 4.054 | 427 |
| 0.03 | full nuclear | 4.640 | 321 | 4.688 | 737 |
| 0.03 | rank-2 truncation | 3.843 | 279 | 3.892 | 737 |

也就是说，按 objective-change 收敛时，真实 F-error 已经略微反弹。

## 对 simulation 的含义

这不意味着正式算法可以用 `M0` 来停止。真实 `M0` 只在 simulation 诊断中可见。

但它说明：

- 只看 penalized objective 的收敛速度，不能完全代表估计误差；
- 需要同时检查 estimator error、rank、singular values 和 objective；
- 对 rank 已知的实验，rank truncation 是一个有意义的后处理诊断。

## 对实际算法的含义

在真实应用中没有 `M0`，因此不能用 true Frobenius error 选择停止点。更现实的做法包括：

- 使用 validation set 或 cross-validation 选择 lambda；
- 检查 estimated rank 和 singular value decay；
- 对多个 lambda 报告稳定区间；
- 在理论允许时使用和 design/noise scaling 匹配的 lambda 公式。
