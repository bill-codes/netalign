%% Evaluate LP solution of grid graphs

%% Setup paths
addpath('../../matlab');
addpath('../misc/gaimc'); % get bfs function
addpath('~/devextern/Clp-1.9.0');

%% Setup parameters
% params(i,:) = [k,q,p]
params = [
    10  1  0.10
    %10  1  0.50
    %10  1  1.00
    20  1  0.025
    20  1  0.05
    %20  1  0.10
    %20  1  0.15
    %20  1  0.20
    %20  2  0.05
    %20  3  0.01
    50  1  0.005
    50  1  0.003
    ];

%%
results = [];
nrep = 10;
alpha = 0; 
beta = 1;
ri = 1;
for i=1:length(params)
    p = params(i,:);
    for ti = 1:nrep
        results(ri).param = p;
        results(ri).k = p(1);
        results(ri).q = p(2);
        results(ri).p = p(3);
        results(ri).rstate = RandStream.getDefaultStream.State;
        [A,B,L,xy] = align_grid_test_data(p(1),p(2),p(3));
        [S,w,li,lj] = netalign_setup(A,B,L);
        [f,C,b,Si,Sj] = netalign_lp_prob(S,w,alpha,beta,li,lj);
        y = clp([],-f,C,b,[],[],zeros(size(f)),ones(size(f)));
        xf = y(1:length(w)); 
        Sf = sparse(Si,Sj,y(length(w)+1:end),size(S,1),size(S,2)); 
        Sf = max(Sf,Sf');
        results(ri).f = f'*y;
        [ma mb mi weight overlap] = mwmround(xf,S,w,li,lj);
        results(ri).round_xf = [overlap weight sum(ma==mb)];
        [ma mb mi weight overlap] = mwmround(S*xf,S,w,li,lj);
        results(ri).round_Sxf = [overlap weight sum(ma==mb)];
        [ma mb mi weight overlap] = mwmround(Sf*xf,S,w,li,lj);
        results(ri).round_Sfxf = [overlap weight sum(ma==mb)];
        
        results(ri)
        ri = ri+1;
    end
end
        
