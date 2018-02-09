% Draw Approximation Errors from different single-pass algorithms.
% Called by CompSinglePass_Plot.m
% It outputs Figure 8 in paper "Efficient randomized algorithms for
%   the fixed-precision low-rank matrix approximation", by W. Yu, et al.

tic
err = [];
k = 200; b = 10;
normA= norm(A1, 'fro');     % for Matrix 1
e1 = SVD_errors(A1, k, b);
err = [err, e1/normA];
[~, ~, ~, e4] = singlePass2011(A1, k, b);
err = [err, e4/normA];
[~, ~, e3] = randQB_FP_k(A1, k, b, 0);
err = [err, e3/normA];
% semilogy(err);

err2=[];
k2 = 100; b = 10;           % for Matrix 2
normA= norm(A2, 'fro');
e1 = SVD_errors(A2, k2, b);
err2 = [err2, e1/normA];
[~, ~, ~, e4] = singlePass2011(A2, k2, b);
err2 = [err2, e4/normA];
[~, ~, e3] = randQB_FP_k(A2, k2, b, 0);
err2 = [err2, e3/normA];
toc

subplot(1,2,1);
semilogy(0:b:k, [1; err(:,1)],  'LineWidth', 1); 
hold on;
xx= (0:b:k)'* ones(1,2);
semilogy((0:b:k), [1 ; err(:,2)], 'x',  ...
    'LineWidth', 1, 'MarkerSize',7);
semilogy((0:b:k), [1 ; err(:,3)], 'ko',   ...
    'LineWidth', 1, 'MarkerSize',7);

legend('SVD', 'single-pass[2]', 'randQB\_FP');        
xlabel('l');
ylabel('||A-QB||/||A||');

subplot(1,2,2);
semilogy(0:b:k2, [1; err2(:,1)],  'LineWidth', 1); 
hold on;
xx= (0:b:k2)'* ones(1,2);
semilogy((0:b:k2), [1 ; err2(:,2)], 'x',  ...
    'LineWidth', 1, 'MarkerSize',7);
semilogy((0:b:k2), [1 ; err2(:,3)], 'ko',   ...
    'LineWidth', 1, 'MarkerSize',7);

legend('SVD', 'single-pass[2]', 'randQB\_FP');        
xlabel('l');
ylabel('||A-QB||/||A||');
