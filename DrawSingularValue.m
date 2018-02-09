% Draw singular values computed by single-pass algorithms.
% Called by CompSinglePass_Plot.m
% It outputs Figure 9 in paper "Efficient randomized algorithms for
%   the fixed-precision low-rank matrix approximation", by W. Yu, et al.

k= 50;
% oversampling parameter is 10 in rSVD_exSP and rSVDbasic and rSVDsp.

subplot(1,2,1);         % for Matrix 1
[~, d1, ~]= singlePass2011_svd(A1, k);
d2 = basicQB_svd(A1, k);
[~, d3, ~] = randQB_FP_svd(A1, k, 10);

semilogy(1:k, d1, 'x', 1:k, d3(1:k), 'ko',...
        1:k, d2, '--', 1:k, s1(1:k), 'b-',  'LineWidth', 1, 'MarkerSize', 4);
legend('single-pass [2]','randQB\_FP', 'randQB','SVD');
(d1-s1(1:k)')./(d3(1:k)-s1(1:k)')
ylabel('\sigma_{j}');

axis([1, 50, 1e-4, 1]);


subplot(1,2,2);         % for Matrix 2
[~, d1, ~]= singlePass2011_svd(A2, k);
d2 = basicQB_svd(A2, k);
[~, d3, ~] = randQB_FP_svd(A2, k, 10);
semilogy(1:k, d1, 'x', 1:k, d3(1:k), 'ko',...
        1:k, d2, '--', 1:k, s2(1:k), 'b-',  'LineWidth', 1, 'MarkerSize', 4);
legend('single-pass [2]','randQB\_FP', 'randQB','SVD');
(d1-s2(1:k)')./(d3(1:k)-s2(1:k)')
ylabel('\sigma_{j}');
    
axis([1, 50, 1e-4, 1]);
