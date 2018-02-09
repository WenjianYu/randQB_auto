function [Q, B, err, errQ, err_id]= basicQB(A, k, P, step)
% [Q, B, err, errQ,Eall]= basicQB(A, k, P, step)
% The basic randQB algorithm. k is the rank parameter
% P is the power parameter, step is an optional argument for outputing 
%    error trend with incresing rank (usually step is factor of k).
% err is the ||A-QB||_F, errQ is ||I-Q'*Q||_inf.
% err_id is the error indicator, ||A||^2-||B||^2.
[m,n]= size(A);
% E= norm(A, 'fro')^2;
E= A(:)'*A(:);    % square of Frobenius norm of A
Omg= randn(n, k);
Y= A*Omg;
[Q, ~]= qr(Y, 0);
for j=1:P,
    [Q, ~]= qr(A'*Q, 0);
    [Q, ~]= qr(A*Q, 0);
end
B= Q'*A;

if nargin>3,          % output error trend with incresing rank
    err=[]; errQ=[]; err_id=[];
    for i=step:step:k,
        err= [err, norm(A-Q(:, 1:i)*B(1:i,:), 'fro')];
        errQ= [errQ, norm(eye(i)- Q(:, 1:i)'*Q(:, 1:i), 'inf')];
%         E= E- norm(B(i-step+1:i,:), 'fro')^2;
        tempB= B(i-step+1:i,:);
        E= E- tempB(:)'*tempB(:);
        err_id= [err_id, E];  
    end
end 