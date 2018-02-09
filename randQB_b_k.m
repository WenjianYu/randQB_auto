function [Q, B, errs] = randQB_b_k(A, k, bs, p)
% [Q, B, errs] = randQB_b_k(A, k, bs, p)
% The randQB_b algorithm with fixed rank k.
% p is the power parameter, bs is rank-increase step (usually a factor of k).
% errs is the error ||A-QB||_F.

    [m, n] = size(A);
    Q = [];
    B = [];
    A_ = A;        % A_ is actually the residual, updated iteratively.
    errs = [];
    i = 1;
    while (i <= k)
        omg = randn(n, bs);
        y = A_ * omg;
        [q, ~] = qr(y, 0);
        for j=1:p,      % added by yuwj
            [q, ~]= qr(A_'*q, 0);
            [q, ~]= qr(A_*q, 0);
        end
        if i > 1
            [q, ~] = qr(q - Q * (Q' * q), 0);
        end
        b = q' * A_;    
        Q = [Q, q];
        B = [B; b];
        A_ = A_ - q * b;
        i = i + bs;
        errs= [errs, norm(A-Q*B, 'fro')];
        %errs = [errs; errors(Q, B, A)];
    end
end