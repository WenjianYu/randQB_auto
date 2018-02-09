function [Q, B, errs, errQ, err_id] = randQB_EI_k(A, k, b, P, opt)
% [Q, B, errs, errQ, Eall] = randQB_EI_k(A, k, b, P, opt)
% The randQB_EI algorithm with fixed rank parameter k.
% P is the power parameter, b is rank-increase step (usually a factor of k).
% opt is a switch for output errors (can be any value).
% errs is the error ||A-QB||_F, errQ is ||I-Q'*Q||_inf.
% err_id is the error indicator, ||A||^2-||B||^2.

    [m, n]  = size(A);
    E= A(:)'*A(:);
    
    Q = zeros(m, 0);
    B = zeros(0, n);
    
    errs = [];
    errQ = [];
    err_id= [];
    r = 1;
    while r < k
        Omg = randn(n, b);
        Y = A * Omg - (Q * (B * Omg));
        [Qi, ~] = qr(Y, 0);
        
        for j = 1:P        % power scheme
            [Qi, ~] = qr(A'*Qi - B'*(Q'*Qi), 0);  % can skip orthonormalization for small b.
            [Qi, ~] = qr(A*Qi - Q*(B*Qi), 0);
        end
        
        if r>1,            % can skip the first re-orthogonalization
            [Qi, ~] = qr(Qi - Q * (Q' * Qi), 0);
        end
        Bi= Qi'*A;  % another choice is Bi = Qi' * A - Qi' * Q * B;
        
        Q = [Q, Qi];
        B = [B; Bi];
        
        r = r + b;
        if nargin>4,    % output errors
            E= E- Bi(:)'*Bi(:);
            errs= [errs, norm(A-Q*B, 'fro')];
            errQ= [errQ, norm(eye(r-1)- Q'*Q, 'inf')];
            err_id= [err_id, E]; 
        end
    end
end