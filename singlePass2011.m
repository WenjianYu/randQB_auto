function [Q, B, Qt, errs] = singlePass2011(A, k, step)
% [Q, B, Qt, errs] = singlePass2011(A, k, step)
% The single-pass algorithm in paper "Finding structure with randomness:
%    Probabilistic algorithms for constructing approximate matrix
%    decompostions" by N. Halko et al.
% k is rank parameter, step is rank-increase step for calculating error trend.
% A is approximated by Q*B*Qt. errs is the approx. error in Frobenius norm.

    [m, n]  = size(A);
    
    errs = [];
    if nargin<3
        r= k;
    else
        r = step;
    end
    while r <= k
        Omg = randn(n, r);
        Omgt= randn(m, r);
        Y= A*Omg; Yt= A'*Omgt;
        [Q, ~]= qr(Y, 0); [Qt, ~]=qr(Yt, 0);
        B= (Omgt'*Q)\(Yt'*Qt);
        
        r = r + step;
        errs = [errs; norm(A- Q*B*Qt', 'fro')];
    end
end