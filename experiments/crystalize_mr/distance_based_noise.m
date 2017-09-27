function N = distance_based_noise(G,q)
% Add edges to G with probability q/(d(u,v)) for all u, v in G
% This implementation is not smart and tries all pairs of edges.

n = size(G,1);
E = zeros(0,2);
for i=1:n
    d = bfs(G,i);         
    p = q./d.^2;
    e = rand(size(p))<=p;   % edges to add
    e(1:i) = 0;             % remove one half of the edges
    newedges = find(e);
    E = [E; [i*ones(length(newedges),1) newedges(:)]];
end
N = sparse(E(:,1),E(:,2),true,size(G,1),size(G,2));
N = N|N';
