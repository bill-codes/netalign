%% Evaluate matching Rameau to LCSH
% Rameau is the french version of LCSH.  These two networks have
% been aligned by hand and we have the matching between them.  Here
% we see how much of that matching automatic translation reproduces.

%% Setup
addpath('../../matlab');

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

%% Output basic facts about the matching
% For each language, how much of the match do we capture there?
fprintf('%35s %6i\n','total matching',nnz(M));
fprintf('%35s %6i %5.1f%%\n','matching captured by L',nnz(L&M),...
    100*nnz(L&M)/nnz(M));
fprintf('%35s %6i %5.1f%%\n','matching captured by L (english)',...
    nnz(Len&M),100*nnz(Len&M)/nnz(M));
fprintf('%35s %6i %5.1f%%\n','matching captured by L (french)',...
    nnz(Lfr&M),100*nnz(Lfr&M)/nnz(M));

%% Network Alignment

%% Run MWM
% THE ALIGNMENT IS NOT A MATCHING!  ARG!
[val ma mb] = bipartite_matching(M);
fprintf('Converting to a matching...\n');
fprintf('Original alignment: %i edges\n',nnz(M));
fprintf('Matching alignment: %i edges (loss %i)\n', ...
    length(ma), nnz(M)-length(ma));
evaluate_alignment(A,B,ma,mb);
Malg = sparse(ma,mb,1,nA,nB);
fprintf('%25s  %i (%5.2f%%)\n', 'correct matches', ...
    nnz(Malg.*M), 100*nnz(Malg.*M)/nnz(M));

%% Try SCBP and BP and MR
objs = [1,1; 1,2; 0,1];
results = [];
results.objs = objs;
for oi=1:length(objs)
    alpha = objs(oi,1);
    beta = objs(oi,2);
    
    % baseline
    t0 = tic;
    mwm = mwmround(w,S,w,li,lj);
    dt = toc(t0);
    results.mwm(oi).time = dt;
    
    curm = w;
    results.mwm(oi).x = curm;
    fprintf('\n\n%25s  (%i,%i)\n','MWM results',alpha,beta);
    [ma mb mi weight overlap] = mwmround(curm,S,w,li,lj);
    results.mwm(oi).weight = weight;
    results.mwm(oi).overlap = overlap;
    [eresults Mcc] = evaluate_alignment(A,B,ma,mb); eresults
    results.mwm(oi).eval = eresults;
    results.mwm(oi).Mcc = Mcc;
    Mbp = sparse(ma,mb,1,nA,nB);
    fprintf('%25s  %i (%5.2f%%)\n', 'correct matches', nnz(Mbp.*M), 100*nnz(Mbp.*M)/nnz(M));
    results.mwm(oi).ncorr = nnz(Mbp.*M);
    results.mwm(oi).ncorr_cc = nnz(Mcc.*M);
    
    
    t0 = tic;
    mbp=netalignbp(S,w,alpha,beta,li,lj,0.999,3,500);
    dt = toc(t0);
    results.bp(oi).x = mbp;
    results.bp(oi).time = dt;
    
    save 'partial.mat' results
    
    curm = mbp;
    fprintf('\n\n%25s  (%i,%i)\n','BP results',alpha,beta);
    [ma mb mi weight overlap] = mwmround(curm,S,w,li,lj);
    results.bp(oi).weight = weight;
    results.bp(oi).overlap = overlap;
    [eresults Mcc] = evaluate_alignment(A,B,ma,mb); eresults
    results.bp(oi).eval = eresults;
    results.bp(oi).Mcc = Mcc;
    Mbp = sparse(ma,mb,1,nA,nB);
    fprintf('%25s  %i (%5.2f%%)\n', 'correct matches', nnz(Mbp.*M), 100*nnz(Mbp.*M)/nnz(M));
    results.bp(oi).ncorr = nnz(Mbp.*M);
    results.bp(oi).ncorr_cc = nnz(Mcc.*M);
    
    t0 = tic;
    mbp=netalignscbp(S,w,alpha,beta,li,lj,0.999,3,500);
    dt = toc(t0);
    results.scbp(oi).x = mbp;
    results.scbp(oi).time = dt;
    
    save 'partial.mat' results
    
    curm = mbp;
    fprintf('\n\n%25s  (%i,%i)\n','SCBP results',alpha,beta);
    [ma mb mi weight overlap] = mwmround(curm,S,w,li,lj);
    results.scbp(oi).weight = weight;
    results.scbp(oi).overlap = overlap;
    [eresults Mcc] = evaluate_alignment(A,B,ma,mb); eresults
    results.scbp(oi).eval = eresults;
    results.scbp(oi).Mcc = Mcc;
    Mbp = sparse(ma,mb,1,nA,nB);
    fprintf('%25s  %i (%5.2f%%)\n', 'correct matches', nnz(Mbp.*M), 100*nnz(Mbp.*M)/nnz(M));
    results.scbp(oi).ncorr = nnz(Mbp.*M);
    results.scbp(oi).ncorr_cc = nnz(Mcc.*M);

    t0 = tic;
    mmr=netalignmr(S,w,alpha,beta,li,lj,0.4,5,[],500);
    dt = toc(t0);
    results.mr(oi).x = mmr;
    results.mr(oi).time = dt;
    
    save 'partial.mat' results
    
    curm = mmr;
    fprintf('\n\n%25s  (%i,%i)\n','MR results',alpha,beta);
    [ma mb mi weight overlap] = mwmround(curm,S,w,li,lj);
    results.mr(oi).weight = weight;
    results.mr(oi).overlap = overlap;
    [eresults Mcc] = evaluate_alignment(A,B,ma,mb); eresults
    results.mr(oi).eval = eresults;
    results.mr(oi).Mcc = Mcc;
    Mbp = sparse(ma,mb,1,nA,nB);
    fprintf('%25s  %i (%5.2f%%)\n', 'correct matches', nnz(Mbp.*M), 100*nnz(Mbp.*M)/nnz(M));
    results.mr(oi).ncorr = nnz(Mbp.*M);
    results.mr(oi).ncorr_cc = nnz(Mcc.*M);

end

save 'rameau2lcsh_results.mat' results
!rm partial.mat


