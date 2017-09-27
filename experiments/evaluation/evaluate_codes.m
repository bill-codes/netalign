%% Evaluate performance of our methods on a problem
% Evaluate a large set of methods with a range of parameters.

%% Setup

experiment = 'evaluation';
[curpath curdir] = fileparts(pwd);
if strcmp(curdir,experiment)==0
    error('experiment:wrongDir',...
        'experiments must be run from the correct directory, should be %s', ...
        experiment);
end        

addpath('../../matlab');

%% Load data
% This section selects which data file we process
if ~exist('prob','var')
    error('must set prob before running this script');
end
[S,w,li,lj,A,B,L] = load_netalign_problem(prob);

%% Methods
% We will evaluate the following methods
%
% isorank with x rounding and a*w + b*S*x rounding
% bp with damping 1 and damping 2
% (for small problems)
%   lp with clp
%

% setup the variables driving the methods
objs = [10,1; 2,1; 1,1; 1,2; 1,10];
alphas = [0.3, 0.5, 0.85, 0.95];
gammas = [0.9 0.99 0.995 0.999];

gammasmr = [0.1 0.4];
stepms = [5 25 50];

%% MWM
mwmresults = [];
[ma mb mi weight overlap] = mwmround(w,S,w,li,lj);
mwmresults(1).hist = [weight length(ma) overlap];
mwmresults(1).prob = prob;
mwmresults(1).str = 'mwm';

save([prob '-mwm.mat'],'mwmresults');

%% IsoRank
isoresults = [];
for ai=1:length(alphas)
    for rtype=1:2
        % the alpha and beta are ignored in this problem
        alpha = alphas(ai);
        dt=tic;
        [x flag hist] = isorank(S,w,1,1,li,lj,alpha,rtype);
        isoresults(end+1).prob=prob;
        isoresults(end).time = toc(dt);
        isoresults(end).alpha = alpha;
        isoresults(end).rtype = rtype;
        isoresults(end).hist = hist(:,3:5);
        isoresults(end).str = sprintf('isorank-r%i-%4.2f',rtype,alpha);
        save([prob '-isorank.mat'],'isoresults');
    end
end

%% Belief Propagation
bpresults = [];
for oi=1:size(objs,1);
    for gi=1:length(gammas)
        for dtype=2:3
            a = objs(oi,1);
            b = objs(oi,2);
            gamma = gammas(gi);
            
            dt=tic;
            [x hista histb] = netalignbp(S,w,a,b,li,lj,gamma,dtype,500);
            
            bpresults(end+1).prob = prob;
            bpresults(end).time = toc(dt);
            bpresults(end).hist = hista(:,2:4);
            bpresults(end).obj = [a,b];
            bpresults(end).gamma = gamma;
            bpresults(end).dtype = dtype;
            bpresults(end).str = sprintf('bp-a-%i-%i-d%i-%4.2f',a,b,dtype,gamma);
            bpresults(end+1)=bpresults(end);
            bpresults(end).hist = histb(:,2:4);
            bpresults(end).str = sprintf('bp-b-%i-%i-d%i-%4.2f',a,b,dtype,gamma);

            save([prob '-bp.mat'],'bpresults');   
        end
    end
end

%% Belief Propagation with Square Constraints
scbpresults = [];
for oi=1:size(objs,1);
    for gi=1:length(gammas)
        for dtype=2:3
            a = objs(oi,1);
            b = objs(oi,2);
            gamma = gammas(gi);
            
            dt=tic;
            [x hista histb] = netalignscbp(S,w,a,b,li,lj,gamma,dtype,500);
            
            scbpresults(end+1).prob = prob;
            scbpresults(end).time = toc(dt);
            scbpresults(end).hist = hista(:,2:4);
            scbpresults(end).obj = [a,b];
            scbpresults(end).gamma = gamma;
            scbpresults(end).dtype = dtype;
            scbpresults(end).str = sprintf('bp-a-%i-%i-d%i-%4.2f',a,b,dtype,gamma);
            scbpresults(end+1)=scbpresults(end);
            scbpresults(end).hist = histb(:,2:4);
            scbpresults(end).str = sprintf('bp-b-%i-%i-d%i-%4.2f',a,b,dtype,gamma);

            save([prob '-scbp.mat'],'scbpresults');   
        end
    end
end

%% Matching Relaxations
mrresults = [];
for oi=1:size(objs,1);
    for stepmi=1:length(stepms)
        for gi=1:length(gammasmr)
            a = objs(oi,1);
            b = objs(oi,2);
            stepm = stepms(stepmi);
            gamma = gammasmr(gi);
            
            % skip a few
            if stepm==50 && gamma==0.4, continue; end
            if stepm==5 && gamma==0.1,continue; end

            dt=tic;
            [x status hist] = netalignmr(S,w,a,b,li,lj,gamma,stepm,1,500);
            

            mrresults(end+1).prob = prob;
            mrresults(end).time = toc(dt);
            mrresults(end).hist = hist(:,5:7);
            mrresults(end).gamma = gamma;
            mrresults(end).stepm = stepm;
            mrresults(end).obj = [a,b];

            mrresults(end).str = sprintf('mr-%i-%i-d%i-%4.2f',a,b,stepm,gamma);


            save([prob '-mr.mat'],'mrresults');   
        end
    end
end

