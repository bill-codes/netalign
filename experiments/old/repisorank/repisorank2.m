function [xr]=repisorank(S,w,a,b,ei,ej,gamma)
%
v = w;
v = v./csum(v);
mi = v;
%gamma = 0.95;

maxiter = 10;
Sk = S;
for i=1:maxiter
    [xr,flag,hist] = isorank(S,v,a,b,ei,ej,[],10,[],normout(Sk));
    [ma,mb,mi,weight,overlap] = mwmround(xr,S,w,ei,ej);
    fprintf('\n');
    fprintf('** RepIsoRank iter %3i : weight= %8g ; overlap = %8i\n', ...
        i, weight, overlap);
    fprintf('\n');
    
    v = gamma*mi + (1-gamma)*v;
    Sk = diag(sparse(v))*Sk*diag(sparse(v));
end
