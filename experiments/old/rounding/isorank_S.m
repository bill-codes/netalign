function [x flag reshist]=isorank_S(S,w,a,b,ei,ej,tol,maxit,verbose,P)
% ISORANK Solve an overlap matching problem with IsoRank
%
% This implementation of IsoRank uses solves the PageRank system using
% the power method and outputs the matching statistics on every iteration.
%
% x=isorank(S,w,a,b,ei,ej) iteratively solves the IsoRank system 
%   (I - b/(a+b) S'*D^+ - corr)x = (1-b/(a+b)) w 
% where corr is just a correction term that fixes the norm of x at 1.
% This system generates a sequence of iterates that should provide
% a good heuristic to solve the QP
%   max  a*w'*x + b/2*x'*S*x
%   s.t. x is a matching in the bipartite graph [ei,ej].
% To generate a matching from an iterate x, we solve a max weight
% bipartite matching problem in the weighted graph [ei,ej,x].  We round
% the solution at every iteration and return the x that generated the 
% best objective value.
%
% x=isorank(S,w,a,b) does not round the solution at every iteration and
% only returns the converged isorank solution.
%
% Typically, we solve this matching problem at each iteration and report
% the change in x, the value of the object function f, the weight of the
% matching w'*rounded(x), the overlap (rounded(x)'*S*rounded(x)/2),
% and the cardinality of the matching (sum(rounded(x)).
%
% The return value is not the matching, to get the matching call mwmround.
%
% Example:
%   load('../data/example-overlap.mat');
%   a = 3; b = 17;
%   x=isorank(S,w,a,b,li,lj); % best solution at any iterate
%   [weight overlap]=mwmround(x,S,w,li,lj);
%   x=isorank(S,w,a,b); % only the converged isorank solution


% David Gleich and Ying Wang
% Copyright, Stanford University, 2008

% 2008-11-03: Initial version based on previous powerpr.m code.

if ~exist('P','var') || isempty(P), P = normout(S); end
v = w./csum(w); assert(all(v>=0)); 

n=size(P,1); 
if ~exist('a','var') || isempty(a), a=0.5; end
if ~exist('b','var') || isempty(b), b=1; end
if ~exist('ei','var'), ei=[]; end
if ~exist('ej','var'), ej=[]; end
if ~exist('tol','var') || isempty(tol), tol=1e-12; end
if ~exist('maxit','var') || isempty(maxit), maxit=100; end
if ~exist('verbose','var') || isempty(verbose), verbose=1; end
em = max(ei); en = max(ej);
allstats = ~isempty(ei) && ~isempty(ej);
if allstats, rhistsize=5; else rhistsize=1; end
r = b/(a+b);
x=zeros(n,1)+v; delta=2; iter=0; reshist=zeros(maxit,rhistsize); 
xbest=x; fbest=0; fbestiter=0;
if verbose && allstats% print the header
    fprintf('%5s   %4s   %8s   %7s %7s %7s %7s\n', ...
        'best', 'iter', 'pr-delta', 'obj', 'weight', 'card', 'overlap');
elseif verbose, fprintf('%4s   %8s\n','iter', 'delta');
end
while iter<maxit && delta>tol
    y=r*(P'*x); omega = 1-csum(y); y = y + omega*v; x = x-y; 
    delta=norm(x,1); reshist(iter+1)=delta; iter=iter+1; x=y*(1/csum(y));
    if allstats
         X = sparse(ei,ej,x,em,en);
        [val ma mb match] = bipartite_matching(a*w + b*S*x,ei,ej,em,en); 
        val=sum(w(match)); overlap=match'*(S*match)/2; 
        f = a*val + b*overlap;
        if f>fbest, xbest=a*w + b*S*x; fbest=f; fbestiter=iter; itermark='*'; 
        else itermark=' ';
        end
    end
    if verbose && allstats
        fprintf('%5s   %4i   %8.1e   %7g %7g %7i %7i\n', ...
            itermark, iter, delta, f, val, length(ma), overlap);
        reshist(iter,2:end) = [a*val + b*overlap, val, length(ma), overlap];
    elseif verbose, fprintf('%4i    %8.1e\n', iter, delta);
    end
end
flag=delta>tol; reshist=reshist(1:iter,:);
if allstats, x=xbest; end

