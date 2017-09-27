%%
addpath('~/bin/mosek/5/toolbox/r2007a/');

%% 
[f,A,b] = netalign_lp_prob(S,w,0,1,li,lj,'sym');

%%
tic;
[y,fval,flag] = linprog(-f,A,b,[],[],zeros(size(f)),ones(size(f)),[],optimset('MaxIter',1e6)); 
toc;
flag
