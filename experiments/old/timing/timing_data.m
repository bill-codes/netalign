
%% Setup

experiment = 'timing';
[curpath curdir] = fileparts(pwd);
if strcmp(curdir,experiment)==0
    error('experiment:wrongDir',...
        'experiments must be run from the correct directory, should be %s', ...
        experiment);
end        

addpath('../../matlab');
%addpath('/home/dgleich/devextern/Clp-1.9.0');

%% Load data

% load('../../private_data/lcsh2wiki-small');
% prob='lcsh2wiki-small';
% w=lw;

% load('../../data/example-2.mat');
% [S,w,li,lj] = netalign_setup(A,B,L);
% prob = 'example-2';

S = readSMAT('../../private_data/lcsubj2wikipedia-qp-squares.smat');
Ldata = load('../../private_data/lcsubj2wikipedia-qp-squares.edges');
li = Ldata(:,1)+1; 
lj = Ldata(:,2)+1;
w = Ldata(:,3);
prob = 'lcsh2wiki-big';

%%

dts = struct;

niter = 100; 

t0=clock();
netalignbp(S,w,1,1,li,lj,[],[],niter);
dts.bp = etime(clock,t0);

t0=clock();
isorank(S,w,1,1,li,lj,0.99,[],0,niter);
dts.isorank = etime(clock,t0);

t0=clock();
netalign_lagrange(S,w,1,1,li,lj,10,[],niter);
dts.llp = etime(clock,t0);

dts

save(sprintf('%s-timing',prob),'dts')

%%

