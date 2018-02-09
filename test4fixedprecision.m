n= 8000;
M1= genTestMatrix(n, n, 1);
M2= genTestMatrix(n, n, 2);
M3= genTestMatrix(n, n, 3);

% test Adpative Range Finder.
% [Q, B, k]= AdpRangeFinder(M1, 1e-2); k
% [Q, B, k]= AdpRangeFinder(M1, 1e-4); k
% [Q, B, k]= AdpRangeFinder(M2, 1e-4); k
% [Q, B, k]= AdpRangeFinder(M2, 1e-5); k
% [Q, B, k]= AdpRangeFinder(M3, 1e-2); k
% [Q, B, k]= AdpRangeFinder(M3, 1.5e-3, 7950); k

% test randQB_EI
tic; [Q, B, k]= randQB_EI_auto(M1, 1e-2, 10, 1); toc
k
tic; [Q, B, k]= randQB_EI_auto(M1, 1e-4, 10, 1); toc
k
tic; [Q, B, k]= randQB_EI_auto(M2, 1e-4, 10, 1); toc
k
tic; [Q, B, k]= randQB_EI_auto(M2, 1e-5, 10, 1); toc
k
tic; [Q, B, k]= randQB_EI_auto(M3, 1e-2, 10, 1); toc
k
tic; [Q, B, k]= randQB_EI_auto(M3, 1.5e-3, 40, 1); toc
k

% test randQB_FP
tic; [Q, B, k]= randQB_FP_auto(M1, 1e-2, 10, 1); toc
k
tic; [Q, B, k]= randQB_FP_auto(M1, 1e-4, 10, 1); toc
k
tic; [Q, B, k]= randQB_FP_auto(M2, 1e-4, 10, 1); toc
k
tic; [Q, B, k]= randQB_FP_auto(M2, 1e-5, 10, 1); toc
k
tic; [Q, B, k]= randQB_FP_auto(M3, 1e-2, 10, 1); toc
k
tic; [Q, B, k]= randQB_FP_auto(M3, 1.5e-3, 40, 1); toc
k


clear all;
readImage;
[Q, B, k]= AdpRangeFinder(A, 0.1); k
tic; [Q, B, k]= randQB_EI_auto(A, 0.1, 10, 1); toc
k
tic; [Q, B, k]= randQB_EI_auto(A, 0.1, 20, 1); toc
k
tic; [Q, B, k]= randQB_EI_auto(A, 0.1, 10, 2); toc
k
tic; [Q, B, k]= randQB_EI_auto(A, 0.1, 20, 2); toc
k
tic; [Q, B, k]= randQB_FP_auto(A, 0.1, 10, 1); toc
k
tic; [Q, B, k]= randQB_FP_auto(A, 0.1, 20, 1); toc
k
tic; [Q, B, k]= randQB_FP_auto(A, 0.1, 10, 2); toc
k
tic; [Q, B, k]= randQB_FP_auto(A, 0.1, 20, 2); toc
k


clear all;
loadAminerMatrix;
% [Q, B, k]= AdpRangeFinder(A, 0.5, 8100); k
tic; [Q, B, k]= randQB_EI_auto(A, 0.5, 50, 1); toc
k
tic; [Q, B, k]= randQB_EI_auto(A, 0.5, 50, 2); toc
k
tic; [Q, B, k]= randQB_FP_auto(A, 0.5, 50, 1); toc
k
tic; [Q, B, k]= randQB_FP_auto(A, 0.5, 50, 2); toc
k


