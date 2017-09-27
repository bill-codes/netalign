function Xlp = netsimlp(A,B)
L = ones(size(A,1),size(B,1));
[S,w,li,lj] = netalign_setup_dir(A,B,L);
[f,C,d] = netsim_lp_prob(S,w,0,1,li,lj);
[x z status] = clp([],-f,C,d,[],[],zeros(length(f),1),ones(length(f),1));
xs = x(1:length(w));
Xlp = full(sparse(li,lj,xs));