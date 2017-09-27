%% Format the Rameau2LCSH results

%% Setup
addpath('../../matlab');
load rameau2lcsh_results;

%% Load data
load('../../data/lcsh2rameau/rameau2lcsh');
whos

%%
% A is rameau (Anames are the names of each node)
% B is lcsh (Bnames are the names of each node)
% M is the exact matching
% A2Ben/fr are matchings from rameau to lcsh in english and in french
% B2Aen/fr are matchings from lcsh to rameau in english and in french

A = A|A';
B = B|B';

nA = size(A,1);
nB = size(B,2);

% fix the sizes of A2Ben
A2Ben(nA,nB)=0;
A2Bfr(nA,nB)=0;
B2Aen(nB,nA)=0;
B2Afr(nB,nA)=0;

%% 
% Setup link graphs
Len = max(A2Ben,B2Aen');
Lfr = max(A2Bfr,B2Afr');
L = max(Len,Lfr);

%%
% Setup network alignment
[S,w,li,lj] = netalign_setup(A,B,L);


%% Baseline
[ma mb mi weight overlap] = mwmround(w,S,w,li,lj);
[mwmeval Mcc] = evaluate_alignment(A,B,ma,mb);

%% Results printing
objs = results.objs;
algs = {'bp','scbp','mr','mwm'};
for oi = 1:length(objs)
    ostr = sprintf('$\\alpha=%i, \\beta=%i$',results.objs(oi,1),results.objs(oi,2));
    for ai = 1:length(algs)
        alg = algs{ai};
        r = results.(alg)(oi);
        fprintf('%15s ', ostr); ostr = '';
        fprintf('& %4s ', upper(alg));
        fprintf('& %7.1f ', r.weight);
        fprintf('& %7i ', r.overlap);
        fprintf('& %6.1f ', r.time);
        fprintf('& %7i ', r.ncorr);
        fprintf('& %5.1f\\%% ', (r.ncorr/nnz(M))*100);
        fprintf('& %5.1f\\%% ', (r.ncorr/r.eval.size)*100);
        fprintf('& %7i ', r.eval.triangles);
        %fprintf('& %7i ', r.eval.largest_component);
        %fprintf('& %7i ', r.ncorr_cc);
        fprintf('\\\\ \n');
    end
    fprintf('\\midrule \\\\ \n');
     
end
        
%% Spit out some details for the text and table
fprintf('Number of exact matches nnz(M) %i\n',nnz(M));
fprintf('Size of Rameau %i\n', size(A,1));
fprintf('Size of LCSH %i\n', size(B,1));
fprintf('Number of possible matches nnz(L) %i\n',nnz(L));
fprintf('Number of captured matches nnz(L.*M) %i\n',nnz(L.*M));
fprintf('Weight of matching sum(sum(L.*M)) %f\n', full(sum(sum(L.*M))));

%%
% Compute the overlap of the exact matching, we need to solve a netalign
% problem to get the answer because M isn't quite a matching
[MS, Mw, mli,mlj] = netalign_setup(A,B,M);
x = netalignmr(MS,Mw,0,1,mli,mlj,[],[],[],100);
[mla mla mlindex weight overlap] = mwmround(x,MS,Mw,mli,mlj);



