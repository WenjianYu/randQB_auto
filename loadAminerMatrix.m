tic;
Acoo= load('Aminer100K_matrix.txt');
m= Acoo(size(Acoo,1), 1)+1;
n= max(Acoo(:,2))+1;
A= sparse(Acoo(:,1)+1, Acoo(:,2)+1, Acoo(:,3), m, n);
clear Acoo;
A= A';      % keywork-person matrix
toc