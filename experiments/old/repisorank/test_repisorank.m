%% Overlap statistics from network alignment algorithms
% This experiment tests overlap numbers between Wikipedia and LCSH 
% from our simple network alignment algorithms.

%% Setup
% Make sure we are in the correct directory
experiment = 'repisorank';
[curpath curdir] = fileparts(pwd);
if strcmp(curdir,experiment)==0
    error('experiment:wrongDir',...
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

alpha = 1; 
beta = 19;
[xiso flag histiso] = isorank(S,w,alpha,beta, li, lj, [], 50);

%%
xrepiso = repisorank(S,w,alpha,beta, li, lj, 0.95);

%%
xrepiso2 = repisorank2(S,w,alpha,beta, li, lj, 0.95);
