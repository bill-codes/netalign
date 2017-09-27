%% Investigate solving the LP using Clp

%% Setup and load
load lcsh2wiki-small
addpath('~/devextern/Clp-1.9.0');
addpath('../../matlab');
name = 'lcsh2wiki-small';
w = lw;

%%
% 
alpha = 0;
beta = 1;
lpform = 'tight';

%%
% load example-overlap
% name = 'example';
%% 
load('../../private_data/lcsh2wiki_full.mat');
name = 'lcsh2wiki-big';

%% 
[f,A,b] = netalign_lp_prob(S,w,alpha,beta,li,lj,lpform);

%%
tic; 
[x,z,status]=clp([],-f,A,b,[],[],zeros(size(f)),ones(size(f)),...
    struct('maxnumseconds',10*86400,'verbose',1));
dt=toc, status
save(sprintf('clp-%s-%s-%i-%i-solution.mat',lpform,name,alpha,beta),'x','z','status','dt'); 

