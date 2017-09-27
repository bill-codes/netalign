%% Evaluate lcsh2wiki small

addpath('~/devextern/Clp-1.9.0');
addpath('~/research/publications/2009/netalign/matlab');
addpath('../rounding/');

load lcsh2wiki-small
w = lw;
nedges = length(w);

format long g;

%%
alpha = 1;
beta = 0;
[ma mb mi mw mo] = mwmround(w,S,w,li,lj); 
fprintf('& %8g ', [alpha*mw + beta*mo, mo, mw, sum(mi)]); fprintf('\\\\\n');
zmwm = mi;

%%
alpha = 1;
beta = 1; 
[f,A,b,Si,Sj] = netalign_lp_prob(S,w,alpha,beta,li,lj,'sym');
[z,y,status]=clp([],-f,A,b,[],[],zeros(size(f)),ones(size(f)));
Sf = sparse(Si,Sj,z(length(w)+1:end),size(S,1),size(S,2));
Sf = max(Sf,Sf');

zlp = z(1:length(w));

zbp = netalignbp(S,w,alpha,beta,li,lj);
zbp_y = netalignbp_y(S,w,alpha,beta,li,lj);
zbp_yz = netalignbp_yz(S,w,alpha,beta,li,lj);

ziso = isorank(S,w,alpha,beta,li,lj,1e-5,50,1);
ziso_S = isorank_S(S,w,alpha,beta,li,lj,1e-5,50,1);

zmwm = mi;
fprintf('& %8g ', f'*z); fprintf('& -------- & -------- & -------- \\\\\n');

[ma mb mi mw mo] = mwmround(zlp,S,w,li,lj); 
fprintf('& %8g ', [alpha*mw + beta*mo, mo, mw, sum(mi)]); fprintf('\\\\\n');

[ma mb mi mw mo] = mwmround(alpha*w + beta/2*full(sum(Sf,2)),S,w,li,lj); 
fprintf('& %8g ', [alpha*mw + beta*mo, mo, mw, sum(mi)]); fprintf('\\\\\n');

[ma mb mi mw mo] = mwmround(alpha*w + beta/2*Sf*zlp,S,w,li,lj); 
fprintf('& %8g ', [alpha*mw + beta*mo, mo, mw, sum(mi)]); fprintf('\\\\\n');

[ma mb mi mw mo] = mwmround(zbp,S,w,li,lj); 
fprintf('& %8g ', [alpha*mw + beta*mo, mo, mw, sum(mi)]); fprintf('\\\\\n');

[ma mb mi mw mo] = mwmround(zbp_y,S,w,li,lj); 
fprintf('& %8g ', [alpha*mw + beta*mo, mo, mw, sum(mi)]); fprintf('\\\\\n');

[ma mb mi mw mo] = mwmround(zbp_yz,S,w,li,lj); 
fprintf('& %8g ', [alpha*mw + beta*mo, mo, mw, sum(mi)]); fprintf('\\\\\n');

[ma mb mi mw mo] = mwmround(ziso,S,w,li,lj); 
fprintf('& %8g ', [alpha*mw + beta*mo, mo, mw, sum(mi)]); fprintf('\\\\\n');

[ma mb mi mw mo] = mwmround(ziso_S,S,w,li,lj); 
fprintf('& %8g ', [alpha*mw + beta*mo, mo, mw, sum(mi)]); fprintf('\\\\\n');

[ma mb mi mw mo] = mwmround(zmwm,S,w,li,lj); 
fprintf('& %8g ', [alpha*mw + beta*mo, mo, mw, sum(mi)]); fprintf('\\\\\n');

%%
alpha = 0;
beta = 1; 
[f,A,b] = netalign_lp_prob(S,w,alpha,beta,li,lj,'sym');
[z,y,status]=clp([],-f,A,b,[],[],zeros(size(f)),ones(size(f)));
zlp = z(1:length(w));

zbp = netalignbp(S,w,alpha,beta,li,lj);
zbp_y = netalignbp_y(S,w,alpha,beta,li,lj);
zbp_yz = netalignbp_yz(S,w,alpha,beta,li,lj);

ziso = isorank(S,w,alpha,beta,li,lj,1e-5,50,1);
ziso_S = isorank_S(S,w,alpha,beta,li,lj,1e-5,50,1);

zmwm = mi;
fprintf('& %8g ', f'*z); fprintf('& -------- & -------- & -------- \\\\\n');

[ma mb mi mw mo] = mwmround(zlp,S,w,li,lj); 
fprintf('& %8g ', [alpha*mw + beta*mo, mo, mw, sum(mi)]); fprintf('\\\\\n');

[ma mb mi mw mo] = mwmround(alpha*w + beta/2*full(sum(Sf,2)),S,w,li,lj); 
fprintf('& %8g ', [alpha*mw + beta*mo, mo, mw, sum(mi)]); fprintf('\\\\\n');

[ma mb mi mw mo] = mwmround(alpha*w + beta/2*Sf*zlp,S,w,li,lj); 
fprintf('& %8g ', [alpha*mw + beta*mo, mo, mw, sum(mi)]); fprintf('\\\\\n');

[ma mb mi mw mo] = mwmround(zbp,S,w,li,lj); 
fprintf('& %8g ', [alpha*mw + beta*mo, mo, mw, sum(mi)]); fprintf('\\\\\n');

[ma mb mi mw mo] = mwmround(zbp_y,S,w,li,lj); 
fprintf('& %8g ', [alpha*mw + beta*mo, mo, mw, sum(mi)]); fprintf('\\\\\n');

[ma mb mi mw mo] = mwmround(zbp_yz,S,w,li,lj); 
fprintf('& %8g ', [alpha*mw + beta*mo, mo, mw, sum(mi)]); fprintf('\\\\\n');

[ma mb mi mw mo] = mwmround(ziso,S,w,li,lj); 
fprintf('& %8g ', [alpha*mw + beta*mo, mo, mw, sum(mi)]); fprintf('\\\\\n');

[ma mb mi mw mo] = mwmround(ziso_S,S,w,li,lj); 
fprintf('& %8g ', [alpha*mw + beta*mo, mo, mw, sum(mi)]); fprintf('\\\\\n');

[ma mb mi mw mo] = mwmround(zmwm,S,w,li,lj); 
fprintf('& %8g ', [alpha*mw + beta*mo, mo, mw, sum(mi)]); fprintf('\\\\\n');
