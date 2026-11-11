## 0. Conventions

### Math Conventions
- A **hat symbol** above a variable means it is an estimated value. **Example:** $`\hat{\Omega}`$ is an estimate for $`\Omega`$
- **Returns are log returns:** $`r_{i,t} = \log X_{i,t} - \log X_{i,t-1}`$. Simple returns $`(X_t - X_{t-1})/X_{t-1}`$ are used only in the robustness check (Section 6). *Reason:* every Section-5 formula is derived in log prices, so known-answer tests are exact.

### Index Conventions
- Asset 0 is cash (`X[:, 0] == 1`). `R[:, 0] == 0` always. Risky assets are `1 ... N-1`.
- Time index `i` runs `0 ... n-1`; return `R[i]` is from price `X[i]` to `X[i+1]`

## 1. Data Layout
| Object | Python Variable | Definition | Shape/type |
|---|---|---|---|
| $`N`$ | `N` | Number of assets **including cash** | int |
| - | `n` | Total number of steps | int |
| $`X`$ | `X` | Price matrix; `X[:, 0]` for cash | `(n+1, N)` float |
| $`r`$ | `R` | Return matrix | `(n, N)` float |
| $`d`$ | `d` | Number of independent Brownian motions | int |
| $`k`$ | `k` | Number of calm portfolios (the true dimension of $`\mathcal{N}`$) | int |
| $`J^A`$ | `J` | Calm Portfolios | `(N, k)` float |
| $`\mu`$ | `mu` | Mean Return Vector (not to be confused with the index $`\mu`$) | `(N,)` float |
| $`\sigma^a_{\mu}`$ | `sigma` | Volatility Matrix | `(N, d)` float |
| $`\alpha_{\mu}`$ | `alpha` | Drift of asset $`\mu`$ | `(N,)` float |
| $`\alpha^A`$ | `alpha_A` | Drift of the null portfolio $`J^A`$ (V&F measure of arbitrage) | `(k,)` float |
| $`\Gamma_0`$ | `Gamma0` | Contemporaneous covariance ($`\Gamma_0[i,j]=\mathrm{Cov}(r_{i,t},r_{j,t})`$) | `(N, N)` float |
| $`\Gamma_1`$ | `Gamma1` | Lag-1 autocovariance ($`\Gamma_1[i,j]=\mathrm{Cov}(r_{i,t+1},r_{j,t})`$) | `(N, N)` float |

## 2. `Series` class
| Attribute | Definition | Shape/type |
|---|---|---|
| `T` | Total amount of time | float |
| `n` | Total number of steps | int |
|  `dt` | Time between steps | float |
| `t` | Array of times | `(n+1,)` float |
