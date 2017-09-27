%% Use the network alignment LP to upper bound the wordnet problem
%
% TODO Describe the wordnet problem

%% Setup
addpath('/home/dgleich/devextern/Clp-1.9.0');
addpath('../../matlab');

%% 
% Load data
L = readSMAT('Ying2/L.mat');
A = readSMAT('Ying2/graph.EN.smat');
B = readSMAT('Ying2/graph.ES.smat');

%%
% Build network alignment problem
[S,w,li,lj]=netalign_setup(A,B,L);

%%
wcount = w;
wones = ones(size(w));

%% 
% Use only overlap
[f,A,b] = netalign_lp_prob(S,w,0,1,li,lj,'sym');

%% Run LP
tic; 
[x,z,status]=clp([],-f,A,b,[],[],zeros(size(f)),ones(size(f)),...
    struct('maxnumseconds',10*86400,'verbose',1));
toc

save 'wordnet-0-1-lp.mat' x z f A b S w li lj
%%
[overlap, matching] = mwmround(x,S,w,li,lj);