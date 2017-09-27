%% Evaluate algorithm performance on a power law graph
% We generate a random power law graph, and evaluate various network
% alignment algorithms on those graphs.

%%
addpath('../misc/gaimc'); % get bfs function
addpath('~/devextern/Clp-1.9.0/');
addpath('../../matlab');

%% compile the generation code
!g++ -O2 RandomPowerLaw.cpp


%% Set the parameters for the experiment

nrep = 48;
theta = 1.8;
n = 400;
q = 0.02;
d = 1;
pns = [linspace(0.5,8,13) 10 12.5 15 20];

alpha = 0;
beta = 1;

%% Construct an area to save data

fvals = zeros(nrep,length(pns));
fvals_ref = zeros(nrep, length(pns));
rvals = zeros(nrep,length(pns));
ncorr = zeros(nrep,length(pns));

rvals_iso = zeros(nrep,length(pns));
ncorr_iso = zeros(nrep,length(pns));

rvals_bp = zeros(nrep,length(pns));
ncorr_bp = zeros(nrep,length(pns));

rvals_scbp = zeros(nrep,length(pns));
ncorr_scbp = zeros(nrep,length(pns));

partialname = sprintf('partial.mat');

matlabpool(3);

for pi = 1:length(pns)
    p = pns(pi)./n;
    % generate the graphs
    cmd = sprintf('./a.out %i %f %f %f %i',n,theta,q,p,nrep);
    system(cmd);
    parfor ri = 1:nrep
        % Generate the test problem
        A = readSMAT(sprintf('A%i.smat',ri-1));
        B = readSMAT(sprintf('B%i.smat',ri-1));
        Ldata = load(sprintf('L%i.data',ri-1));
        Ldata(:,1)=Ldata(:,1)+1; Ldata(:,2)=Ldata(:,2)+1;
        Ldata(end+1,:)=[n,n,0];
        L = spconvert(Ldata);
        A = A|A';
        B = B|B';
        
        % Setup the problem
        [S,w,li,lj] = netalign_setup(A,B,L);
        
        % Compute the reference solution
        [ma mb mi weight overlap] = mwmround(li==lj,S,w,li,lj);
        fvals_ref(ri,pi) = alpha*weight + beta*overlap;
        
        % Solve with the matching relaxation
        [xmr status] = netalignmr(S,w,alpha,beta,li,lj,0.4,25,[],5000);
        [ma mb mi weight overlap] = mwmround(xmr,S,w,li,lj);
        rvals(ri,pi) = alpha*weight + beta*overlap;
        ncorr(ri,pi) = sum(ma==mb);        
        fvals(ri,pi) = status(3); % get the upper bound on the solution
        rv = rvals(ri,pi); % save the rounded objective and number of correct matches
        nc = ncorr(ri,pi);
        
        % compute the solution with belief propagation
        xbp = netalignbp(S,w,1,2,li,lj,0.999,[],500);
        [ma mb mi weight overlap] = mwmround(xbp,S,w,li,lj);
        rvals_bp(ri,pi) = alpha*weight + beta*overlap;
        ncorr_bp(ri,pi) = sum(ma==mb);
        
        % compute solutions with enhanced belief propagation
        xbp = netalignscbp(S,w,1,2,li,lj,0.999,[],500);
        [ma mb mi weight overlap] = mwmround(xbp,S,w,li,lj);
        rvals_scbp(ri,pi) = alpha*weight + beta*overlap;
        ncorr_scbp(ri,pi) = sum(ma==mb);
        
        % compute with Isorank
        xiso = isorank(S,w,1,19,li,lj);
        [ma mb mi weight overlap] = mwmround(xiso,S,w,li,lj);
        rvals_iso(ri,pi) = alpha*weight + beta*overlap;
        ncorr_iso(ri,pi) = sum(ma==mb);
    end
    save(partialname, 'pi','pns','fvals','fvals_ref','rvals','ncorr','*_scbp','*_bp','*_iso')
end
save(sprintf('results-n-%i-nrep-%i.mat',n,nrep), 'pns','fvals','fvals_ref','rvals','ncorr','*_scbp','*_bp','*_iso')
delete(partialname)

!rm A*.smat 
!rm B*.smat 
!rm S*.smat 
!rm L*.data


