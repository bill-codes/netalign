%% Investigate solving the LP using SNOPT

%% Setup and load
load lcsh2wiki-small
addpath('~/devextern/snopt7/matlab');
w = lw;

%% 
[f,A,b] = netalign_lp_prob(S,w,1,1,li,lj,'sym');
nvar = length(f);
ncon = length(b);
Flow = [-Inf; -Inf*ones(ncon,1)];
Fupp = [Inf; b];
xlow = zeros(nvar,1);
xupp = ones(nvar,1);
ObjAdd = 0;
ObjRow = 1;
Ahat = [-f'; A];
[iAfun jAvar Aval] = find(Ahat);

%%
x0 = zeros(nvar,1);

%%
snopt_netalign_lpfun([],ncon);

%%
snset ('Defaults');
snprint ('netalign-lp-solution.out');
snspec ('netalign-lp.spc');
[x,F,inform,xmul,Fmul]=snopt(x0,xlow,xupp,Flow,Fupp,'snopt_netalign_lpfun',...
                        ObjAdd,ObjRow,Aval,iAfun,jAvar,[],[]);
                    
%%

