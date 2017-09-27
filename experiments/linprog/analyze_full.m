%% Run some quick analyses on the results from the full matching problem

load_big_lcsh2wiki
load clp-lcsh2wiki-0-1-full-solution
y = x;
x = x(1:length(w));
[f,C,b] = netalign_lp_prob(S,w,0,1,li,lj);
[f,C,b,Si,Sj] = netalign_lp_prob(S,w,0,1,li,lj);
Sf = sparse(Si,Sj,y(length(w)+1:end),size(S,1),size(S,2));
Sf = max(Sf,Sf');

%% Value histograms

set(0,'defaultaxesfontsize',12,'defaultaxeslinewidth',.7,...
'defaultlinelinewidth',.8,'defaultpatchlinewidth',.7) 

semilogy(sort(x),'k.-','MarkerSize',6);
set(gca,'XTick',[]); ylabel('value'); grid on;
print('-depsc2','lcsh2wiki-full-0-1-xvals.eps');


semilogy(sort(y(length(w)+1:end)),'k.-','MarkerSize',6);
set(gca,'XTick',[]); ylabel('value'); grid on;
print('-depsc2','lcsh2wiki-full-0-1-yvals.eps');



%%
format long g;
[ma mb mi weight overlap] = mwmround(x, S, w, li, lj); 
[overlap weight sum(mi)]
%%
[ma mb mi weight overlap] = mwmround(S*x, S, w, li, lj);
[overlap weight sum(mi)]
%%
[ma mb mi weight overlap] = mwmround(full(sum(Sf,2)), S, w, li, lj);
[overlap weight sum(mi)]
%%
[ma mb mi weight overlap] = mwmround(Sf*x, S, w, li, lj);
[overlap weight sum(mi)]
%% 
Sf2 = spfun(@(x) x.*(x>=1e-5),Sf);
x2 = full(spfun(@(x) x.*(x>=1e-5),x));
%%
[overlap weight] = mwmround(x2, S, w, li, lj)
%%
[overlap weight] = mwmround(S*x2, S, w, li, lj)
%%
[overlap weight] = mwmround(full(sum(Sf2,2)), S, w, li, lj)
%%
[overlap weight] = mwmround(Sf2*x2, S, w, li, lj)

%% Reanalyze the other methods

load_big_lcsh2wiki
load clp-lcsh2wiki-0-1-full-solution
z = x;
[f,C,b] = netalign_lp_prob(S,w,0,1,li,lj);
[f,C,b,Si,Sj] = netalign_lp_prob(S,w,0,1,li,lj);
Sf = sparse(Si,Sj,z(length(w)+1:end),size(S,1),size(S,2));
Sf = max(Sf,Sf');

zlp = z(1:length(w));

zbp = netalignbp(S,w,alpha,beta,li,lj,0.85,50);
zbp_y = netalignbp_y(S,w,alpha,beta,li,lj,0.85,50);
zbp_yz = netalignbp_yz(S,w,alpha,beta,li,lj,0.85,50);

ziso = isorank(S,w,alpha,beta,li,lj,1e-5,50,1);
ziso_S = isorank(S,w,alpha,beta,li,lj,1e-5,50,1);

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
