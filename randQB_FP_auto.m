function [Q, B, k] = randQB_FP_auto(A, relerr, b, P)
% [Q, B, k] = randQB_FP_auto(A, relerr, b, P)
% The fixed-precision randQB_FP algorithm.
% It produces QB factorization of A, whose approximation error fulfills
%     ||A-QB||_F <= ||A||_F* relerr.
% b is block size, P is power parameter.
% Output k is the rank.

    [m, n]  = size(A);
    Q = zeros(m, 0);
    B = zeros(0, n);

    maxiter= 50;                % this may be changed case by case.
    maxiter= min(maxiter, ceil(min(m,n)/3/b));
    l= b*maxiter;               % an emperical setting of l.
    Omg = randn(n, l);
    while P > 0     % power scheme
        [G, ~] = qr(A * Omg, 0);
        [Omg, ~] = qr(A' * G, 0);
        P = P - 1;
    end
    G = A * Omg;
    H = A' * G;
    % =========
    E= norm(A, 'fro')^2;
    E0= E;
    threshold= relerr^2*E;    
    r = 1;
    flag= false;
    
    for i=1:maxiter,
        t = B * Omg(:, r:r+b-1);
        Y = G(:, r:r+b-1) - (Q * t);
        [Qi, R] = qr(Y, 0);
        
        [Qi, R1] = qr(Qi - Q * (Q' * Qi), 0);
        R = R1 * R;
           
        Bi = R' \ ((H(:, r:r+b-1))' - (Y' * Q) * B- t'*B);
        
        Q = [Q, Qi];
        B = [B; Bi];
        r = r + b;
        
        temp = E- norm(Bi, 'fro')^2;
        
        if temp< threshold,     % for precise rank determination 
            for j=1:b,
                E= E-norm(Bi(j,:))^2;
                if E< threshold,
                    flag= true;
                    break;
                end
            end
        else
            E= temp;
        end
        if flag,
            k= (i-1)*b+j;
            break;
        end
    end
    if ~flag,
        fprintf('E = %f. Fail to converge within maxiter!\n', sqrt(E/E0));
    end
end
