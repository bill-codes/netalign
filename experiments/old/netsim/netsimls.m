function [X,M,rhs]=netsimls(A,B,L,a,b)
m = size(A,1);
n = size(B,1);
if ~exist('L','var') || isempty(L), L = ones(m,n); end;
if ~exist('a','var') || isempty(a), a = 0; end;
if ~exist('b','var') || isempty(b), b = 1; end;
M = [b/2*kron(B',A')                 -kron(ones(n,1),speye(m));
     kron(ones(1,n),speye(m))       sparse(m,m)];
rhs = [-a*L(:); ones(m,1)];
x = M\rhs;
X = reshape(x(1:m*n),m,n);
