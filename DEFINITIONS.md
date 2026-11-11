## 1. Data Layout
| Object | Definition | Shape/type |
|---|---|---|
| $`N`$ | Number of assets **including cash** | int |

**Returns are log returns:** $`r_{i,t} = \log X_{i,t} - \log X_{i,t-1}`$. Simple returns $(X_t - X_{t-1})/X_{t-1}$ are used only in the robustness check (Section 6). *Reason:* every Section-5 formula is derived in log prices, so known-answer tests are exact.

## 2. `Series` class
| Attribute | Definition | Shape/type |
|---|---|---|
| `T` | Total amount of time | float |
| `n` | Total number of steps | int |
|  `dt` | Time between steps | float |
| `t` | Array of times | `(n+1,)` float |

### Example Code:
