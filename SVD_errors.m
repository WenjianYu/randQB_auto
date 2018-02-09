function [errs, Q, B] = SVD_errors(A, k, b)
% [errs, Q, B] = SVD_errors(A, k, b)
% Using standard svd to get the optimal approximation error.
% k is the rank parameter, b is rank-increase step (usually a factor of k).
% errs is the approximation error in Frobenius norm.

    [U, S, V] = svd(A);
    errs = [];
    for i = 1:k/b
        Q = U(:, 1:i*b);
        B = S(1:i*b, 1:i*b) * V(:, 1:i*b)';      
        errs = [errs; norm(A-Q*B, 'fro')];
    end