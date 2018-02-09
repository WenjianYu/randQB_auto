% EIplot.m for validationg Error indicator in randQB_EI.
% It produces the Figure 2 in paper "Efficient randomized algorithms for
%   the fixed-precision low-rank matrix approximation", by W. Yu, et al.
%% for subfigure(a)
[A, d]= gen_rand_mat_exp_decay(3000, 3000, 20);
k=500;
[Q, B, err1, errQ1,sqB1]= basicQB(A, k, 0, 10);
[Q, B, err4, errQ4, sqB4] = randQB_EI_k(A, k, 10, 0, 1);
normA2= norm(A, 'fro')^2;
semilogy(10:10:k, sqB1/normA2, 'b.', 10:10:k, sqB4/normA2, 'ro');
hold on;
semilogy(10:10:k, err1.^2/normA2, 'b-', 'LineWidth', 1);
semilogy(10:10:k, err4.^2/normA2, 'r-', 'LineWidth', 1);
semilogy(10:10:k, errQ1, 'b--', 'LineWidth', 1);
semilogy(10:10:k, errQ4, 'r--', 'LineWidth', 1);
axis([150, 450, 1e-18, 1e-5]);
legend('error indicator(randQB)','error indicator(randQB\_EI)', 'error^2 (randQB)', ...
'error^2(randQB\_EI)', '||Q*Q-I|| (randQB)', '||Q*Q-I|| (randQB\_EI)');
xlabel('l');
% ylabel('normalized square of error / error indicator');
h= gcf;
set(h,  'Position',[560 100 400 420]);

%% for subfigure (b) 
[A, d]= gen_rand_mat_exp_decay(20000, 20000, 200);
k=4000;

[Q, B, err1, errQ1,sqB1]=basicQB(A, k, 0, 10);
[Q, B, err4, errQ4, sqB4] = randQB_EI_k(A, k, 10, 0, 1);
normA2= norm(A, 'fro')^2;
semilogy(10:10:k, sqB1/normA2, 'b.', 10:10:k, sqB4/normA2, 'ro');
hold on;
semilogy(10:10:k, err1.^2/normA2, 'b-', 'LineWidth', 1);
semilogy(10:10:k, err4.^2/normA2, 'r-', 'LineWidth', 1);
semilogy(10:10:k, errQ1, 'b--', 'LineWidth', 1);
semilogy(10:10:k, errQ4, 'r--', 'LineWidth', 1);
axis([3200, 4000, 1e-17, 1e-11]);
legend('error indicator(randQB)','error indicator(randQB\_EI)', 'error^2 (randQB)', ...
'error^2(randQB\_EI)', '||Q*Q-I|| (randQB)', '||Q*Q-I|| (randQB\_EI)');
xlabel('l');
% ylabel('normalized square of error / error indicator');
h= gcf;
set(h,  'Position',[560 100 400 420]);