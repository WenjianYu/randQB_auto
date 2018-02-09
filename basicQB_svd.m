function [U, S, V]=basicQB_svd(A, k, P)
% [U, S, V]=basicQB_svd(A, k, P)
% Rank-k truncated SVD of A based on the basic randQB algorithm.
% Syntax:
%  s = basicQB_svd(A, k)
%  s = basicQB_svd(A, k, P)
%  [U, S, V]= basicQB_svd(A, k)
%  [U, S, V]= basicQB_svd(A, k, P)
%  -P is an optional parameter to balance time and accuacy (default value 0).
%   With large P, the accuracy increases with runtime overhead.

if nargin<3,
    P=0;
end

s=10;               % over-sampling
[m,n]= size(A);
B= randn(n, k+s);
U= A*B;
[U, ~]= qr(U, 0);
for j=1:P,
    [B, ~]= qr(A'*U, 0);   % May reduce an orthogonalization
    [U, ~]= qr(A*B, 0);
end
B= A'*U;

if nargout==1,
    U= svd(B','econ');
    U= U(1:k);
else
    [U1, S, V]= svd(B', 'econ');
    U= U*U1(:,1:k);
    S= diag(S);
    S= S(1:k);
    V= V(:,1:k);
end

end 