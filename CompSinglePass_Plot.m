% Comparing different single-pass algorithms.
% Produce Figure 8 and 9 in paper "Efficient randomized algorithms for
%   the fixed-precision low-rank matrix approximation", by W. Yu, et al.


n= 2000;
[A1, s1] = genTestMatrix(n, n, 1);
[A2, s2] = genTestMatrix(n, n, 2);

DrawApproxError;

figure(2);
DrawSingularValue;
