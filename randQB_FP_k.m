function [Q, B, errs] = randQB_FP_k(A, k, bsize, p)
% [Q, B, errs] = rankQB_FP_k(A, k, bsize, p)
% The randQB_FP algorithm with fixed rank k.
% p is the power parameter, bsize is block size (usually a factor of k).
% errs is the error ||A-QB||_F.
    [m, n]  = size(A);
    Q = zeros(m, 0);
    B = zeros(0, n);
    
    Omg = randn(n, k);
    while p > 0
        [G, ~] = qr(A * Omg, 0);
        [Omg, ~] = qr(A' * G, 0);
        p = p - 1;
    end
    G = A * Omg;
    H = A' * G;
    % =========
    errs = [];
    r = 1;
    while r < k
        t = B * Omg(:, r:r+bsize-1);
        y = G(:, r:r+bsize-1) - (Q * t);
        [q, R] = qr(y, 0);
        
        [q, R1] = qr(q - Q * (Q' * q), 0);
        R = R1 * R;
        b = R' \ ((H(:, r:r+bsize-1))' - (y' * Q) * B- t'*B);
        
        Q = [Q, q];
        B = [B; b];
        r = r + bsize;
        errs = [errs; norm(A-Q*B, 'fro')];
    end
end
