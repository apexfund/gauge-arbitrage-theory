import numpy as np

# Orthonormal basis of zero-cost portfolios over (asset 1, asset 2, cash)
E = np.column_stack([np.array([1, -1, 0]) / np.sqrt(2),
                     np.array([1,  1, -2]) / np.sqrt(6)])

def predict(G0, G1):
    G = np.zeros((3, 3)); G[:2, :2] = G0          # cash has zero variance
    vals, vecs = np.linalg.eigh(E.T @ G @ E)      # eigenvalues ascending
    J = E @ vecs[:, 0]                            # lowest-variance unit mix
    var = vals[0]
    D = J[:2] @ G1 @ J[:2]                        # cash has zero return
    return {"J": J, "var": var, "D_pred": D, "rho_pred": D / var}