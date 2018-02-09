# randQB_auto
randomized QB factorization for fixed-precision low-rank matrix approximation.

This package includes Matlab codes for the randQB_EI and randQB_FP algorithms.
They are efficient randomized algorithms for the fixed-precision low-rank matrix 
approximation. The test cases and scripts for running the experiments in paper
"Efficient randomized algorithms for the fixed-precision low-rank matrix approximation" 
by Wenjian Yu, Yu Gu and Yaohang Li, are also included.


1. Main algorithms

randQB_EI_auto.m -- fixed-precision version of the randQB_EI algorithm

randQB_FP_auto.m -- fixed-precision version of the randQB_EI algorithm

randQB_EI_k.m -- fixed-rank version of the randQB_EI algorithm

randQB_FP_k.m -- fixed-rank version of the randQB_EI algorithm

randQB_FP_svd.m -- compute rank-k truncated SVD with the randQB_FP algorithm

2. Auxiliary algorithms for comparison

basicQB.m -- the basic randQB algorithm (fixed-rank) in [1]

randQB_b_k.m -- the blocked randQB algorithm (fixed-rank) in [2]

AdpRangeFinder.m -- adaptive randomized range finder algorithm (fixed-precision) [1]

singlePass2011.m -- the single-pass algorithm in [1]

singlePass2011_svd.m -- compute rank-k truncated SVD with the single-pass algorithm in [1]

basicQB_svd.m -- compute rank-k truncated SVD with the basic randQB algorithm [1]

SVD_errors.m -- compute the optimal rank-k approximation error with SVD.

3. Test data and codes

genTestMatrix.m -- generate the three dense test matrices (Matrix 1/2/3).

gen_rand_mat_exp_decay.m -- generate a matrix with singular value decay exponentially.

image1.jpg -- A scenic image

Aminer100K_matrix.txt -- A keyword-person matrix from "AMiner" (in COO format)

Aminer100K_s.mat -- Accurate singular values of Aminer100K matrix (obtained with SVD)

4. Experiment scripts.

EIplot.m -- For validating the error indicator in randQB_EI. Also draw Fig. 2 in [3].

CompSinglePass_Plot.m -- Compare different single-pass algorithms. Draw Fig. 8/9 in [3].

DrawSinglarValue.m -- Needed by CompSinglePass_Plot.m

DrawApproxError.m -- Needed by CompSinglePass_Plot.m

test4fixedprecision.m -- Validate the algorithms for fixed-precision computation.

readImage.m -- Convert an image to a matrix.

loadAminerMatrix.m -- Load the Aminer matrix.


For any comments, please contact Wenjian Yu (yu-wj@tsinghua.edu.cn). Thanks!

(Sorry, due to its large size Aminer100K_matrix.txt can't be uploaded here. If need it, please contact me)

Reference

[1] N. Halko, P.-G. Martinsson and J. A. Tropp, "Finding structure with randomness: Probabilistic algorithms for constructing approximate matrix decompositions," SIAM Review, 53 (2011), no. 2, pp. 217{288.

[2] P.-G. Martinsson and S. Voronin, "A randomized blocked algorithm for efficiently computing rank-revealing factorizations of matrices," SIAM J. Sci. Comput., 38(2016), no. 5, pp. S485 - S507.

[3] Wenjian Yu, Yu Gu, and Yaohang Li, "Efficient randomized algorithms for the fixed-precision low-rank matrix approximation," submitted to SIMAX.
