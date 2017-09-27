%% Generate the example matrices used in the paper
load('../../data/example-overlap');
% permute nodes to get labels in the paper
[S,w,li,lj] = netalign_setup(A,B,L);
[f,C,b] = netalign_lp_prob(S,w,0,1,li,lj);
Cmatch = C(1:(size(A,1)+size(B,1)),1:length(li));
lA = [2 1 3 4 5 6];
lB = [2 1 3 4 5];

%% 
% The matrix A
full(Cmatch)

%%
% The matrix S
full(S)

%%
% The edge list
for i=1:length(li)
    fprintf('(%i,%i'')\n',lA(li(i)),lB(lj(i)));
end

%% Solutions
% MWM
[ma mb mi weight overlap] = mwmround(w,S,w,li,lj);
weight
overlap

for i=1:length(ma)
    fprintf('%i <-> %i''\n',lA(ma(i)),lB(mb(i)));
end

%%
% In this case, it matches 1' to 6 and 2' to 1, 

%%
% netalign
x = netalignbp(S,w,0,1,li,lj);
[ma mb mi weight overlap] = mwmround(x,S,w,li,lj);


for i=1:length(ma)
    fprintf('%i <-> %i''\n',lA(ma(i)),lB(mb(i)));
end
