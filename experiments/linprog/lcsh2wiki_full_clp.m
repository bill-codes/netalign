%% Investigate solving the LP using Clp

%% Setup and load
load lcsh2wiki-small
addpath('~/devextern/Clp-1.9.0');
addpath('../../matlab');
w = lw;
%%
%load example-overlap
%% 
load_big_lcsh2wiki

%% 
[f,A,b] = netalign_lp_prob(S,w,1,1,li,lj,'sym');

%%
tic; 
[x,z,status]=clp([],-f,A,b,[],[],zeros(size(f)),ones(size(f)),...
    struct('maxnumseconds',10*86400,'verbose',1));
dt=toc, status
save 'clp-lcsh2wiki-1-1-full-solution.mat' x z status dt; 

