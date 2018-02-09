function [Q, B, k]= AdpRangeFinder(A, relerr, maxcol)
% [Q, B, k]= AdpRangeFinder(A, relerr)
% The adaptive randomized range finder algorithm in paper "Finding structure 
%    structure with randomness: Probabilistic algorithms for constructing 
%    approximate matrix decompostions" by N. Halko et al.
% It produces QB factorization of A, whose approximation error fulfills
%     ||A-QB||_F <= ||A||_F* relerr.
% maxcol is optional, meaning the maxium number of columns tested.
% Parameter r is set to 10, as suggested in that paper.
% Output k is the rank.

[m, n]=size(A);
if nargin<3,
    maxcol= 5000;    % for precalculating a large number of random vectors.
end
tic;
Omg= randn(n, maxcol);
Ymax= A*Omg;
time1= toc;

tic;
r= 10;
y= Ymax(:,1:r);
normy= zeros(r, 1);
for i=1:r,
    normy(i)= norm(y(:,i));
end

j=0;
Q=zeros(m, 0);
acc= relerr*norm(A, 'fro')/(10*sqrt(2/pi));
while max(normy(j+1:j+r))> acc && r+j< maxcol,
    j= j+1;
    y(:,j)= y(:,j)-Q*(Q'*y(:,j));
    q= y(:,j)/norm(y(:,j));
    Q= [Q q];
    newy= Ymax(:, r+j);
    newy= newy- Q*(Q'*newy);
    y=[y newy];
    normy= [normy; norm(newy)];
    y(:,j+1:j+r-1)= y(:,j+1:j+r-1)- q*(q'*y(:,j+1:j+r-1));
end
if max(normy(j+1:j+r))> acc
    fprintf('Warning: the accuracy is not attained with the maximum iteration: %d', maxcol);
end

B= Q'*A;
time2=toc;
k= size(Q, 2);

totalTime= time2+ time1/maxcol*(k+r)

