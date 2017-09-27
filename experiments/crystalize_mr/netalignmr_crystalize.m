%% Test netalignmr crystalization
% We observed that netalignmr often makes progress with the objective value
% when the step size is nearly machine precision.  In this test, we add a
% forced `crystalization' phase to the netalignmr code that always reduces
% the step size to machine precision before returning.  The hope is that
% this procedure might improve its performance on the grid_test examples.
%
% Results:
%   It makes netalignmr a bit better, but not enough to be worth it.  

%
k = 20;
q = 2;
d = 1;
p = 15/k^2; % use an expected degree of 15 -- a nasty case for netalignmr.

%% Construct data
[A,B,L] = align_grid_test_data(k,q,p,d);
L = spones(L); % ignore weights for this test
[S,w,li,lj] = netalign_setup(A,B,L);

alpha=0;
beta=1;

% compute solutions with the matching relaxation
[xmr status] = netalignmr(S,w,alpha,beta,li,lj,0.4,25,[],1000);
% compute the solution 
[xmr2 status2] = netalignmr2(S,w,alpha,beta,li,lj,0.4,25,[],1000);
%