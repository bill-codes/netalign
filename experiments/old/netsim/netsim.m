%% Metrics for network similarity (not alignment)
% Basically, we just drop the matching constraint on the columns of L so
% that only one node must be similar to another one.

%% Setup paths
addpath('~/devextern/Clp-1.9.0/');
addpath('../../matlab');
%% Simple example
%
B = [0 1; 0 0];
A = [0 1 1 0
     1 0 0 1
     0 0 0 1
     0 0 1 0];
Xlp = netsimlp(A,B); 
%%
% Compare against Van Dooren's algorithm
Xi = netsimiter(A,B)

%% 
% Solve with a linear system

%% Try an undirected graph
B = [0 1; 0 0];
A = [0 1 1 0
     1 0 0 1
     0 0 0 1
     0 0 1 0];
A = A|A';

L = ones(size(A,1),size(B,1));
[S,w,li,lj] = netalign_setup_dir(A,B,L);
[f,C,d] = netsim_lp_prob(S,w,0,1,li,lj);
[x z status] = clp([],-f,C,d,[],[],zeros(length(f),1),ones(length(f),1));
xs = x(1:length(w));
Xlp = full(sparse(li,lj,xs))
%%
% Compare against Van Dooren's algorithm
Xi = netsimiter(A,B)     

%%