function X = netsimiter(A,B,tol)
if ~exist('tol','var') || isempty(tol), tol=1e-4; end;
X = ones(size(A,1),size(B,1));
for i=1:100
    Y = A*X*B' + A'*X*B;
    Y = Y./norm(Y,inf);
    delta = norm(Y-X);
    X = Y;
    if delta<tol,
        break;
    end
end
if delta>tol,
    warning('netsimiter:didNotConverge',...
        'the iteration failed to converge to tolerance %s',delta);
end