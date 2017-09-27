%% Linear Programs for Network Alignment
% In this experiment, we investigate MILP formulations of the Network 
% Alignment problem and their solutions.
%
% 2009-02-27: This experiment is done, we have the LP formulation.  The
% real script is netalign_lp_prob

load example-overlap


%% Solve exactly 
xe = netalign_exact(S,w,1,1,li,lj)

%% Binary optimization formulation
nedges = length(li); 
m = max(li); 
n = max(lj);

nsquares = nnz(S);
[si sj] = find(S);
S1 = sparse(si,1:nsquares,1,nedges,nsquares);
S2 = sparse(sj,1:nsquares,1,nedges,nsquares);

f = [w(:); ones(nsquares,1)];

Ar = sparse(li,1:nedges,1,m,nedges); 
Ac = sparse(lj,1:nedges,1,n,nedges); 

Amatch = [Ar sparse(m,nsquares); Ac sparse(n,nsquares)];
bmatch = ones(m+n,1);

Aval = [-S1' speye(nsquares); -S2' speye(nsquares)];
bval = zeros(2*nsquares,1);

A = [Amatch; Aval];
b = [bmatch; bval];

[z,fval,flag] = bintprog(-f,A,b,[],[],[],optimset('MaxIter',1e6)); flag

%% Symmetric binary optimization formulation
% This formulation difference from the one above it because the matrix $Y$
% is constrainted to be symmetic and additional constraints eliminated.
% The problem setup is more complicated, but it has half the number of
% variables.
%
% Variables: nedges + nsquares
% Linear Constraints: n+m + 2*nsquares
% Bound 

nedges = length(li); 
m = max(li); 
n = max(lj);

nsquares = nnz(S)/2;
[si sj] = find(triu(S));
S1 = sparse(si,1:nsquares,1,nedges,nsquares);
S2 = sparse(sj,1:nsquares,1,nedges,nsquares);

f = [w(:); 2*ones(nsquares,1)];

Ar = sparse(li,1:nedges,1,m,nedges); 
Ac = sparse(lj,1:nedges,1,n,nedges); 

Amatch = [Ar sparse(m,nsquares); Ac sparse(n,nsquares)];
bmatch = ones(m+n,1);

Aval = [-S1' speye(nsquares); -S2' speye(nsquares)];
bval = zeros(2*nsquares,1);

A = [Amatch; Aval];
b = [bmatch; bval];

[z,fval,flag] = bintprog(-f,A,b,[],[],[],optimset('MaxIter',1e6)); flag

%% RelaxedLP
nedges = length(li); 
m = max(li); 
n = max(lj);

nsquares = nnz(S);
[si sj] = find(S);
S1 = sparse(si,1:nsquares,1,nedges,nsquares);
S2 = sparse(sj,1:nsquares,1,nedges,nsquares);

f = [w(:); ones(nsquares,1)];

Ar = sparse(li,1:nedges,1,m,nedges); 
Ac = sparse(lj,1:nedges,1,n,nedges); 

Amatch = [Ar sparse(m,nsquares); Ac sparse(n,nsquares)];
bmatch = ones(m+n,1);

Aval = [-S1' speye(nsquares); -S2' speye(nsquares)];
bval = zeros(2*nsquares,1);

A = [Amatch; Aval];
b = [bmatch; bval];

[y,fval,flag] = linprog(-f,A,b,[],[],0,1,[],optimset('MaxIter',1e6)); flag


%% Symmetric Relaxed LP 
% Now we'll look at relaxing the integer constraint.

nedges = length(li); 
m = max(li); 
n = max(lj);

nsquares = nnz(S)/2;
[si sj] = find(triu(S));
S1 = sparse(si,1:nsquares,1,nedges,nsquares);
S2 = sparse(sj,1:nsquares,1,nedges,nsquares);

f = [w(:); 2*ones(nsquares,1)];

Ar = sparse(li,1:nedges,1,m,nedges); 
Ac = sparse(lj,1:nedges,1,n,nedges); 

Amatch = [Ar sparse(m,nsquares); Ac sparse(n,nsquares)];
bmatch = ones(m+n,1);

Aval = [-S1' speye(nsquares); -S2' speye(nsquares)];
bval = zeros(2*nsquares,1);

A = [Amatch; Aval];
b = [bmatch; bval];

[y,fval,flag] = linprog(-f,A,b,[],[],zeros(size(f)),ones(size(f)),[],optimset('MaxIter',1e6)); flag


