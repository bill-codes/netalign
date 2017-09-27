%% Evaluate algorithm performance on a grid using only the SCBP algorithm
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
%
% In this experiment, we only evaluate the SCBP algorithm to try to
% determine if changing any of the rounding modes improves its performance.
% 

%%
addpath('../misc/gaimc'); % get bfs function
addpath('~/devextern/Clp-1.9.0/');
addpath('../../matlab');


%% Set the parameters for the experiment

nrep = 24;

k = 20;
q = 2;
d = 1;
% use expected degree to indicate noise in L
pns2 = [8 10 12.5 15];


%% Run the experiment
alpha = 0;
beta = 1;

fvals_ref2 = zeros(nrep, length(pns2));

rvals_scbp2 = zeros(nrep,length(pns2));
ncorr_scbp2 = zeros(nrep,length(pns2));

for pi = 1:length(pns2)
    p = pns2(pi)./k^2;
    for ri = 1:nrep
        [A,B,L] = align_grid_test_data(k,q,p,d);
        L = spones(L); % ignore weights for this test
        [S,w,li,lj] = netalign_setup(A,B,L);
        
        % compute the reference values
        [ma mb mi weight overlap] = mwmround(li==lj,S,w,li,lj);
        fvals_ref2(ri,pi) = alpha*weight + beta*overlap;
        
        % compute solutions with enhanced belief propagation
        % gamma,dtype,maxiter
        xbp = netalignscbp(S,w,1,2,li,lj,0.999,2,500);
        [ma mb mi weight overlap] = mwmround(xbp,S,w,li,lj);
        rvals_scbp2(ri,pi) = alpha*weight + beta*overlap;
        ncorr_scbp2(ri,pi) = sum(ma==mb);
    end
end

%%

load 'results-k-20-nrep-48.mat'

clf;
hs=plot(...,
    pns, mean(ncorr_bp)/400, 'k.-', ...
    pns, mean(ncorr_scbp)/400, 'b.-', ...
    pns2, mean(ncorr_scbp2)/400, 'r.-');
xlim([0,20]); ylim([0,1.1]);
line(xlim,[1 1],'LineWidth',0.5,'Color',0.5*[1,1,1]);
set(hs,'MarkerSize',13);
ylabel('fraction correct');
hl = legend('BP','SCBP','SCBP2',...
    'Location','SW');
legend('boxoff');
xlabel('expected degree of noise in L (p \cdot n)');
