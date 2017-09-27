function [xr]=repisorank(S,w,a,b,ei,ej,gamma)
%
v = w;
v = v./csum(v);
mi = v;
%gamma = 0.95;

maxiter = 10;

for i=1:maxiter
    v = gamma*mi + (1-gamma)*v;
    [xr,flag,hist] = isorank(S,v,a,b,ei,ej,[],30);
    [ma,mb,mi,weight,overlap] = mwmround(xr,S,w,ei,ej);
    fprintf('\n');
    fprintf('** RepIsoRank iter %3i : weight= %8g ; overlap = %8i\n', ...
        iter, weight, overlap);
    fprintf('\n');
    
    v = gamma*mi + (1-gamma)*v;
end

an
