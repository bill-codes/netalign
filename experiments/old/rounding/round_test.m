%% Analyze rounding modes on small grids

addpath('~/devextern/Clp-1.9.0');
addpath('~/dev/matlab-bgl-4.0');
addpath('../grid_test'); % we need grid test functions

%% Setup a problem that is slightly difficult for the LP solver
% We'll call align_grid_test_data with these parameters
k = 10;
q = 1;
p = 0.1;
alpha = 0;
beta = 1;


%%
nrep = 10;
results = cell(nrep,1);
parfor ri=1:nrep
    result = struct();
    [A,B,L] = align_grid_test_data(k,q,p);
    [S,w,li,lj] = netalign_setup(A,B,L);
    [f,C,b,Si,Sj] = netalign_lp_prob(S,w,alpha,beta,li,lj);
    y = clp([],-f,C,b,[],[],zeros(size(f)),ones(size(f)));
    xf = y(1:length(w)); 
    Sf = sparse(Si,Sj,y(length(w)+1:end),size(S,1),size(S,2)); 
    Sf = max(Sf,Sf');
    fval = f'*y;
    result.f = fval;
        
    [ma mb mi weight overlap] = mwmround(xf,S,w,li,lj);
    rv = (alpha*weight + beta*overlap)/fval;
    nc = sum(ma==mb)/100;
    result.lp = [rv nc];
    
    [ma mb mi weight overlap] = mwmround(S*xf,S,w,li,lj);
    rv = (alpha*weight + beta*overlap)/fval;
    nc = sum(ma==mb)/100;
    result.lp_y = [rv nc];
    
    [ma mb mi weight overlap] = mwmround(Sf*xf,S,w,li,lj);
    rv = (alpha*weight + beta*overlap)/fval;
    nc = sum(ma==mb)/100;
    result.lp_zy = [rv nc];
    
    r = isorank(S,w,alpha,beta,li,lj,1e-5,50,0);
    [ma mb mi weight overlap] = mwmround(r,S,w,li,lj);
    rv = (alpha*weight + beta*overlap)/fval;
    nc = sum(ma==mb)/100;
    result.isorank = [rv nc];
    
    r = isorank_S(S,w,alpha,beta,li,lj,1e-5,50,0);
    [ma mb mi weight overlap] = mwmround(r,S,w,li,lj);
    rv = (alpha*weight + beta*overlap)/fval;
    nc = sum(ma==mb)/100;
    result.isorank_S = [rv nc];
    
    mbest = netalignbp(S,w,alpha,beta,li,lj,0.85,50,0);
    [ma mb mi weight overlap] = mwmround(mbest,S,w,li,lj);
    rv = (alpha*weight + beta*overlap)/fval;
    nc = sum(ma==mb)/100;
    result.bp = [rv nc];
    
    mbest = netalignbp_y(S,w,alpha,beta,li,lj,0.85,50,0);
    [ma mb mi weight overlap] = mwmround(mbest,S,w,li,lj);
    rv = (alpha*weight + beta*overlap)/fval;
    nc = sum(ma==mb)/100;
    result.bp_y = [rv nc];
    
    mbest = netalignbp_yz(S,w,alpha,beta,li,lj,0.85,50,0);
    [ma mb mi weight overlap] = mwmround(mbest,S,w,li,lj);
    rv = (alpha*weight + beta*overlap)/fval;
    nc = sum(ma==mb)/100;
    result.bp_yz = [rv nc];
    
    result
    
    results{ri} = result;
    
end

%% process results

clear sresults;
for ri=1:nrep
    sresults(ri) = results{ri};
end

%%
mresults = [];
fn = fieldnames(sresults(1));
for fi =1:length(fn)
    mresults.(fn{fi}) = mean(...
        reshape([sresults.(fn{fi})], length(sresults(1).(fn{fi})), nrep),2)';
end
% print table
for fi =2:length(fn)
    fprintf('\\prog{%10s} & %5.2f\\%% & %5.2f\\%% \\\\ \n', fn{fi}, 100*mresults.(fn{fi}));
end
%%

    