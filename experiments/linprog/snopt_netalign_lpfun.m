function [F,G] = snopt_netalign_lpfun(x,ncon_in)

persistent ncon z
if isempty(x)
    ncon = ncon_in;
    z = zeros(ncon,1);
    F = [];
    G = [];
    return
end

F = [0; z];
G = [];
