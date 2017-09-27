load('../../private_data/lcsh2wiki_full.mat');
addpath('../../matlab');

[x,status,hist] = netalignmr(S,w,0,1,li,lj,0.4,25,[],10000)
save 'lcsh2wiki_full_mr.mat' x status hist;
