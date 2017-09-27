%% Evaluate algorithm performance on a grid
% In this experiment, we evaluate the performance of the matching
% algorithms on two corrupted grids.
%
% We use
%   [A,B,L] = align_grid_test_data(k,q,p,d) 
% to construct two k-by-k grid graph corrupted by noise in the form of
% edges sampled with probability q/d(u,v)^2 where d(u,v) is the number of
% links between u and v.  These two graphs are A and B.
% To construct L, we first fix the true alignment between A and B as the
% identity between the underlying grids.  Next, if d>0, we randomly add
% edges from all vertices within distance d of the exact match, 
% samping them with probability proporational to their distance.  
% (In all of these tests, d=1.)  Finally, we add random matches to L with
% probability p.

%%
addpath('../misc/gaimc'); % get bfs function
addpath('~/devextern/Clp-1.9.0/');
addpath('../../matlab');


%% Set the parameters for the experiment

nrep = 48;

k = 20;
q = 2;
d = 1;
% use expected degree to indicate noise in L
pns = [linspace(0.5,8,13) 10 12.5 15 20];


%% Run the experiment
alpha = 0;
beta = 1;
fvals = zeros(nrep,length(pns));
rvals = zeros(nrep,length(pns));
ncorr = zeros(nrep,length(pns));

fvals_ref = zeros(nrep, length(pns));

rvals_iso = zeros(nrep,length(pns));
ncorr_iso = zeros(nrep,length(pns));

rvals_bp = zeros(nrep,length(pns));
ncorr_bp = zeros(nrep,length(pns));

rvals_scbp = zeros(nrep,length(pns));
ncorr_scbp = zeros(nrep,length(pns));


matlabpool(3);

for pi = 1:length(pns)
    p = pns(pi)./k^2;
    parfor ri = 1:nrep
        [A,B,L] = align_grid_test_data(k,q,p,d);
        L = spones(L); % ignore weights for this test
        [S,w,li,lj] = netalign_setup(A,B,L);
        
        % compute the reference values
        [ma mb mi weight overlap] = mwmround(li==lj,S,w,li,lj);
        fvals_ref(ri,pi) = alpha*weight + beta*overlap;
        
        % compute solutions with the matching relaxation
        [xmr status] = netalignmr(S,w,alpha,beta,li,lj,0.4,25,[],5000);
        
        % round the output
        [ma mb mi weight overlap] = mwmround(xmr,S,w,li,lj);
        
        % compute the rounded function values
        rvals(ri,pi) = alpha*weight + beta*overlap; 
        ncorr(ri,pi) = sum(ma==mb);        
        fvals(ri,pi) = status(3);
        rv = rvals(ri,pi);
        nc = ncorr(ri,pi);

        % compute solutions with belief propagation
        xbp = netalignbp(S,w,1,2,li,lj,0.999,[],1000);
        [ma mb mi weight overlap] = mwmround(xbp,S,w,li,lj);
        rvals_bp(ri,pi) = alpha*weight + beta*overlap;
        ncorr_bp(ri,pi) = sum(ma==mb);
        
        % compute solutions with enhanced belief propagation
        xbp = netalignscbp(S,w,1,2,li,lj,0.999,[],1000);
        [ma mb mi weight overlap] = mwmround(xbp,S,w,li,lj);
        rvals_scbp(ri,pi) = alpha*weight + beta*overlap;
        ncorr_scbp(ri,pi) = sum(ma==mb);
        
        % compute solutions with 
        xiso = isorank(S,w,1,19,li,lj);
        [ma mb mi weight overlap] = mwmround(xiso,S,w,li,lj);
        rvals_iso(ri,pi) = alpha*weight + beta*overlap;
        ncorr_iso(ri,pi) = sum(ma==mb);

        disp([pi pns(pi) ri status(3) rv nc ...
            rvals_bp(ri,pi) ncorr_bp(ri,pi) ...
            rvals_iso(ri,pi) ncorr_iso(ri,pi)])
    end
    
    save(sprintf('results-k-%i-nrep-%i-partial.mat',k,nrep), 'pi','pns','fvals','fvals_ref','rvals','ncorr','*_scbp','*_bp','*_iso')
end
save(sprintf('results-k-%i-nrep-%i.mat',k,nrep), 'pns','fvals','fvals_ref','rvals','ncorr','*_scbp','*_bp','*_iso')
return

