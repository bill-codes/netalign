
dir = '~yw1984/Experiments/Lagrangean/natalie-0.9';
Ad = load(fullfile(dir,'Mus_m.data'));
Bd = load(fullfile(dir,'Homo_s.data'));
Ld = load(fullfile(dir,'L.data'));

A = sparse(Ad(:,1)+1,Ad(:,2)+1,1);
nA = max(size(A));
A(nA,nA)=0;
B = sparse(Bd(:,1)+1,Bd(:,2)+1,1);
nB = max(size(B));
B(nB,nB)=0;

L = sparse(Ld(:,1)+1,Ld(:,2)+1,1,nA,nB);

A = A-diag(diag(A));
B = B-diag(diag(B));
A=A|A';
B=B|B';

addpath('~/research/publications/2009/netalign/matlab/');
addpath('~/devextern/Clp-1.9.0/');
[S,w,li,lj] = netalign_setup(A,B,L);


%%
alpha=0;
beta=1;

[f,C,b,Si,Sj] = netalign_lp_prob(S,w,alpha,beta,li,lj);
y = clp([],-f,C,b,[],[],zeros(size(f)),ones(size(f)));