%% Network alignment algs for bioinformatics data
% In this experiment, we test our ideas on the original bioinformatics
% data that Singh et al. used in their IsoRank paper.
%
% *This experiment was largely informative.*  See the
% evaluation/evaluate_codes for the actual comparison.

%% Setup

experiment = 'bioinfo';
[curpath curdir] = fileparts(pwd);
if strcmp(curdir,experiment)==0
    error('experiment:wrongDir',...
        'experiments must be run from the correct directory, should be %s', ...
        experiment);
end        

addpath('../../matlab');
addpath('/home/dgleich/devextern/Clp-1.9.0');

%%
% Load files
A = readSMAT('dmela.smat');
A = max(A,A');
B = readSMAT('scere.smat');
B = max(B,B');
L = readSMAT('dmela-scere.smat');
Liso = readSMAT('dmela-scere-isorank-match.smat'); % read the isorank scores from Singh et al.
[S,w,li,lj] = netalign_setup(A,B,L);


%% Run the algs on the sparse matching

%%
% BP algorithm
xbp = netalignbp(S,w,1,1,li,lj);
%%
% Isorank
xiso = isorank(S,w,0.1,0.9,li,lj);
%%
% Lp
[f,C,b] = netalign_lp_prob(S,w,0,1,li,lj,'sym');
tic; 
[x,z,status]=clp([],-f,C,b,[],[],zeros(size(f)),ones(size(f)),...
    struct('maxnumseconds',10*86400,'verbose',1));
toc

%% Expand the matching and then run it
L2 = expand_match(A,B,L);
L2s= -1*spones(L2) + 2*spones(L);
[S,w,li,lj] = netalign_setup(A,B,L2s);
w(w<0.5) = 0;

%%
% Isorank
xiso = isorank(S,w,0.1,0.9,li,lj);

%%
% BP algorithm
xbp = netalignbp(S,w,1,1,li,lj,0.99);

%%
% BP algorithm with randomly perturbed weights
w2=w;
w2(w==0) = min(w(w>0))*1e-2*rand(sum(w==0),1);
xbp = netalignbp(S,w2,1,1,li,lj,0.99);

%% Expand the again matching and then run it
L3 = expand_match(A,B,L2);
[S,w,li,lj] = netalign_setup(A,B,L3);

%%
% Isorank
xiso = isorank(S,w,0.1,0.9,li,lj);

%% Run the algs on the isorank matching subset
% This set can be a bit tricky to form due to memory issues, so we take
% care here to make everything work smoothly.
%
% Actually, it really doesn't work without significant swapping, so we 
% only consider a top-k subset
[li lj lv] = find(Liso);
[ignore p] = sort(lv,1,'descend');
N = 1e5;
Liso2 = sparse(li(p(1:N)), lj(p(1:N)), lv(p(1:N)), size(Liso,1),size(Liso,2));
[Se Le] = make_squares(A,B,Liso2);

%%
spparms('chmodsp',1);

%%
li = Le(:,1);
lj = Le(:,2);
w = Le(:,3);
clear Le;
Se1 = Se(:,1);
Se2 = Se(:,2);
clear Se;
%%
S = sparse(Se1,Se2,true,length(li),length(li));
clear Se1 Se2;

%%
netalignmbp(S,w,1,1,li,lj);

%% Run the algs on the dense matching
Lfull = -ones(size(L))+spones(L)+L;
[Sf,wf,lfi,lfj] = netalign_setup(A,B,Lfull);
wf(w<0.5) = 0; % reset entries to 0

%%
%
xiso = isorank(Sf,wf,0.1,0.9,lfi,lfj);


