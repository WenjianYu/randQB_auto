function [A, d] = genTestMatrix(m, n, t)
% [A, d] = genTestMatrix(m, n, t)
% Generate Matrix 1/2/3 in paper "Efficient randomized algorithms for
%   the fixed-precision low-rank matrix approximation", by W. Yu, et al.
% A is the mxn matrix. t is type identity, 1, 2, or 3.
    L = randn(m, m);
    [U, ~] = qr(L);
    L = randn(n, n);
    [V, ~] = qr(L);
    p= min(m,n);
    d= zeros(1, p);
    g = rand(1, p);
    
    switch t,
        case 1 
            d= (1:p).^(-2);
        case 2 
            d= exp(-(1:p)/7);
        case 3 
            d= 0.0001+1./(1+exp((1:p)-30));
        otherwise
            disp('Not a valid t value!');
            exit;
    end

    S= spdiags(d', 0, m, n);
    A = U * S * V;

