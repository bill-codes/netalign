load('../../data/example-2.mat')
G = [A L; L' B];
Lbig = spaugment(L,0);
Axy_fig = [Axy(:,1) Axy(:,2)+0.1];
Bxy_fig = [Bxy(:,1)+7 Bxy(:,2)-.25];
ABxy = [Axy_fig; Bxy_fig];
gplot(G,ABxy,'.-'); axis equal;

write_tkz_graph(A,Axy_fig,[],0,0);
write_tkz_graph([sparse(13,13) sparse(13,12); sparse(12,13) B],ABxy,[],0,0);
write_tkz_graph(Lbig,ABxy,[],0,0);