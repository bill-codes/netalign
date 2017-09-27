A_e = load('~/Experiments/Lagrangean/natalie-0.9/Mus_m.data');
B_e = load('~/Experiments/Lagrangean/natalie-0.9/Homo_s.data');
L_e = load('~/Experiments/Lagrangean/natalie-0.9/L.data');

A = sparse(A_e(:,1)+1,A_e(:,2)+1,1);
B = sparse(B_e(:,1)+1,B_e(:,2)+1,1);
L = sparse(L_e(:,1)+1,L_e(:,2)+1,1);

A(max(size(A,1),size(A,2)),max(size(A,1),size(A,2))) = 0;
B(max(size(B,1),size(B,2)),max(size(B,1),size(B,2))) = 0;

A = A|A';
B = B|B';
