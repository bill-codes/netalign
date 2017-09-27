function [A,B,L,xy,G,Lstart,Lrand] = align_grid_test_data(k,q,p,d,nm)
% Produce 
%
% Call:
%   align_grid_test_data(k,q,p) 
%   align_grid_test_data(k,q,p,[d],[nm])
%
% The default value of d=0 and nm=0
%
% d is the m
%
% The output is a text problem
%
% A = k-by-k grid with random edges (u,v) sampled with probability q/(d(u,v))
% B = k-by-k grid with random edges (u,v) sampled with probability q/(d(u,v))
% L = exact match between k-by-k grid with random edges with probability p
%     and biased with distance noise out to d links from the main match
%     and a noisy match score if nm=1

if ~exist('d','var') || isempty(d), d=0; end
if ~exist('nm','var') || isempty(nm), nm=0; end

e = ones(k,1); 
GL = spdiags([e,e],[-1,1],k,k);
GI = speye(k);
G = kron(GL,GI) + kron(GI,GL);
G = G - diag(diag(G));
[X,Y] = meshgrid(-1:2/(k-1):1,-1:2/(k-1):1);
xy = [X(:) Y(:)];
A = G+distance_based_noise(G,q);
B = G+distance_based_noise(G,q);
if nm==1
    % create a noisy match
    d = 1+0.05*randn(k^2,1); % values along diagonal
    d = min(d,1); % at most 1
    d(d<0.5) = 0.5; % min at 0.5;
    L = diag(sparse(d));
else
    L = speye(k^2);
end
Lstart = L;

if d>0
    Ld = speye(k^2);
    Ldall = 0;
    for i=1:d
        Ld = expand_match(A,B,Ld);
        Ldall = Ldall + Ld;
    end    
    % sample from Ld with probabiltiy portional to the distance
    [li lj lp] = find(Ldall);
    maxp = max(lp);
    keepl = rand(length(lp),1)<(lp./maxp);
    Ld = sparse(li(keepl),lj(keepl),lp(keepl)./maxp,k^2,k^2);
    L = max(L,Ld);
end
Lrand = sprand(k^2,k^2,p);
L = L + Lrand;
