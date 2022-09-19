# L1 Tucker Decomposition
This repository contains alternating optimization algorithms for L1-norm tucker decomposition

Dependencies: Tensor Toolbox (https://tensortoolbox.com/)

## Algorithms:

**L1HOSVD**
  - **Usage:** [U, G, Xhat, L1met, stats]=L1HOSVD(X, Ks, Uin, *'tol'*, tol, *'maxit'*, maxit, *'X_clean'*, X_clean, *'Un_true'*, Un_true)
  - **Inputs:**
    - **X:** I-way Input Tensor of size (D = [D1 D2 ... DI])
    - **Ks:** Low-rank tensor core size (Ks = [d1 d2 ... dI])
    - **Uin:** Initial Factor matrices (cell array of length I, where Uin{i} is a matrix of size Di by di)
    - *(Optional)* **tol:** optimization tolerance (default value: 1e-4)
    - *(Optional)* **maxit:** maximum number of algorithm iterations  (default value: 50)
    - *(Optional)* **X_clean:** Nominal input tensor (input tensor without noises or outliers) [If given, stats.RERR returns reconstruction error vs. basis updates]
    - *(Optional)* **Un_true:** Nominal factor matrices (If given, algorithm returns subspace error vs. basis update indexes, in stats.SERR)
  - **Outputs:**
    - **U:** Factor matrices
    - **G:** Tensor core    (tensor of size Ks)
    - **Xhat:** Reconstruction of X
    - **L1met:** The value of the objective function (L1-metric of tensor core) at the end of computation
    - **stats:** struct containing some information about algorithm execution
      - stats.update_types : list of the indexes of basis (0 for updating B), in the order they are updated in algorithm execution
      - stats.RERR: list of reconstruction errors after each basis update
      - stats.SERR: list of subspace errors after each basis update

**L1HOOI function:**
  - **Usage:** [U,G,Xhat,L1met_end,stats,funcname,stats_T1] = L1HOOI(X, Ks, Uin, *'tol'*, tol, *'maxit'*, maxit, *'X_clean'*, X_clean, *'Un_true'*, Un_true, 'T',inf, 'selection', selection)
  - **Inputs:**
    - **X:** I-way Input Tensor of size (D = [D1 D2 ... DI])
    - **Ks:** Low-rank tensor core size (Ks = [d1 d2 ... dI])
    - **Uin:** Initial Factor matrices (cell array of length I, where Uin{i} is a matrix of size Di by di)
    - *(Optional)* **proj:** 'L2'(default) for L1HOOI with L2-projection, 'L1' for L1HOOI with L1-projection
    - *(Optional)* **tol:** optimization tolerance (default value: 1e-8)
    - *(Optional)* **maxit:** maximum number of algorithm iterations  (default value: 1000)
    - *(Optional)* **X_clean:** Nominal input tensor (input tensor without noises or outliers) [If given, stats.RERR returns reconstruction error vs. basis updates]
    - *(Optional)* **Un_true:** Nominal factor matrices (If given, algorithm returns subspace error vs. basis update indexes, in stats.SERR)
    - *(Optional)* **T:** If set to a number, it controls the maximum iterations for the outer loop( e.g. set T = 1 for one round of U matrix updates)
    - *(Optional)* **selection:** input options are 'default' (for increasing order), 'random' (update in a random order of basis matrices randomly picked at each outer iteration), or a list specifying the order of basis updates (e.g. [3, 1, 2])
  - **Outputs:**
    - **U:** Factor matrices
    - **G:** Tensor core    (tensor of size Ks)
    - **Xhat:** Reconstruction of X
    - **L1met_end: The value of the objective function (L1-metric of tensor core) at the end of computation
    - **stats:** struct containing some information about algorithm execution
      - stats.update_types : list of the indexes of basis (0 for updating B), in the order they are updated in algorithm execution
      - stats.RERR: list of reconstruction errors after each basis update
      - stats.SERR: list of subspace errors after each basis update
    - **funcname** name of the function configuration (e.g. 'L1HOOI/L2proj': is returned for above function call if T is not given)
    - **stats_T1:** struct containing some information about algorithm execution, until the end of round 1
    
## L1HOOI Configurations:
  - **L1-HOOI/L2 projection (Increasing Order):** [U,G,Xhat] = L1HOOI(X, Ks, Uin)
  - **L1-HOOI/L1 projection (Increasing Order):** [U,G,Xhat] = L1HOOI(X, Ks, Uin, 'proj', 'L1')
  - **L1-HOOI/L2 projection (random permutation Order):** [U,G,Xhat] = L1HOOI(X, Ks, Uin, 'selection', 'random')
  - **L1-HOOI/L1 projection (random permutation Order):** [U,G,Xhat] = L1HOOI(X, Ks, Uin, 'proj', 'L1', 'selection', 'random')
  - **L1-HOOI/L2 projection (T=1):** [U, G, Xhat] = L1HOOI(X, Ks, Uin, 'T', 1) (or use  [\~ ,\~ ,\~ ,\~ ,\~ ,\~ , stats_T1] = L1HOOI(X, Ks, Uin) for concurrent computation of L1-HOOI/L2 projection and L1-HOOI/L2 projection(T=1))
  - **L1-HOOI/L1 projection (T=1):** [U, G, Xhat] = L1HOOI(X, Ks, Uin, 'proj', 'L1', 'T', 1)



