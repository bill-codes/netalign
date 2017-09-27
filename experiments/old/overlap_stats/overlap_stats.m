%% Overlap statistics from network alignment algorithms
% This experiment tests overlap numbers between Wikipedia and LCSH 
% from our simple network alignment algorithms.

%% Setup
% Make sure we are in the correct directory
experiment = 'overlap_stats';
[curpath curdir] = fileparts(pwd);
if strcmp(curdir,experiment)==0
    error('experment:wrongDir',...
        'experiments must be run from the correct directory, should be %s', ...
        experiment);
end        

datadir = fullfile('..','..','private_data');

%%
% Now load our data sets
try
    S = readSMAT(fullfile(datadir,'lcsubj2wikipedia-qp-squares.smat'));
    Ldata = load(fullfile(datadir,'lcsubj2wikipedia-qp-squares.edges'));
    li = Ldata(:,1)+1; 
    lj = Ldata(:,2)+1;
    w = Ldata(:,3);
catch
    fprintf('%s\nTry adding netalign/matlab to your Matlab path?\n%s\n',...
        repmat('*',1,10),repmat('*',1,10));
    rethrow(lasterr)
end

%%
% Setup the experiment sets
experiments = [];

experiments(end+1).label = '1-19';
experiments(end).alpha = 1;
experiments(end).beta = 19;

experiments(end+1).label = '1-2';
experiments(end).alpha = 1;
experiments(end).beta = 2;

experiments(end+1).label = '1-1';
experiments(end).alpha = 1;
experiments(end).beta = 1;

experiments(end+1).label = '2-1';
experiments(end).alpha = 2;
experiments(end).beta = 1;

experiments(end+1).label = '19-1';
experiments(end).alpha = 19;
experiments(end).beta = 1;

%% Run the experments on BP and IsoRank
figure(1); clf; curlegend = {};
hold all;
for ei=1:length(experiments)
    exp=experiments(ei)
    alpha = exp.alpha; beta = exp.beta; label=exp.label;
    
    [xbp hista histb] = netalignbp(S,w,alpha,beta,li,lj,[],50);
    curlegend{end+1}=sprintf('bp-a-%s',label);
    curlegend{end+1}=sprintf('bp-b-%s',label);
    plot(hista(:,3),hista(:,4),'.'); plot(histb(:,3),histb(:,4),'.'); 
    legend(curlegend{:}); drawnow;
    
    [xiso flag histiso] = isorank(S,w,alpha,beta,li,lj,[],50);
    curlegend{end+1}=sprintf('isorank-%s',label);
    plot(histiso(:,4),histiso(:,5),'.'); 
    legend(curlegend{:}); drawnow;
    
    opts = {'delimiter',' ','precision','%18g'};
    dlmwrite(sprintf('isorank-%s.hist',label),histiso,opts{:});
    dlmwrite(sprintf('bp-a-%s.hist',label),hista,opts{:});
    dlmwrite(sprintf('bp-b-%s.hist',label),histb,opts{:});
end
