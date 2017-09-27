%% Plot grid_test results


%% Setup
addpath('../tools/');

get_colors;

set(0,'defaultaxesfontsize',18,'defaultaxeslinewidth',1.5,...
'defaultlinelinewidth',2,'defaultpatchlinewidth',.7) 


%% Load data

%%
% partial
load 'partial.mat';
outputname = 'partial';

%%
% k=20, nrep=48
n=400; nrep=48;
load(sprintf('results-n-%i-nrep-%i.mat',n,nrep));
outputname=sprintf('powerlaw-n-%i-nrep-%i',n,nrep);

%%

clf
hs = plot(...
    pns, mean(fvals./fvals_ref), 'k-', ...     % upper mr
    pns, mean(rvals./fvals_ref), 'k.-', ...    % rounded mr
    pns, mean(rvals_bp./fvals_ref), 'b.-', ...  % rounded bp
    pns, mean(rvals_scbp./fvals_ref), 'r.-', ...  % rounded bp
    pns, mean(rvals_iso./fvals_ref), 'g.-' ...% rounded iso
    );
xlim([0,20]); ylim([0,1.1]);    
line(xlim,[1 1],'LineWidth',0.5,'Color',0.5*[1,1,1]);
set(hs,'MarkerSize',13);
set(hs(3),'Color',niceblue);
set(hs(4),'Color',nicered);
set(hs(5),'Color',nicegreen);
ylabel('rounded objective values');
hl = legend('MR-upper','MR','BP','BPSC','IsoRank',...
    'Location','SW');
legend('boxoff');
% legend('lp solution', 'rounded lp solution', 'correct matches lp', ...
%     'correct matches bp',...
%     'Location','NW');
xlabel('expected degree of noise in L (p \cdot n)');

greygrid;
set(gca,'Color','none'), set(gcf,'Color','none');
set(gcf, 'InvertHardCopy','off');
print(gcf,'-depsc2','-cmyk',[outputname '-fval.eps']);

%%
clf;
hs=plot(...,
    pns, mean(ncorr)/400, 'k.-', ...
    pns, mean(ncorr_bp)/400, 'b.-', ...
    pns, mean(ncorr_scbp)/400, 'r.-', ...
    pns, mean(ncorr_iso)/400, 'g.-');
xlim([0,20]); ylim([0,1.1]);
line(xlim,[1 1],'LineWidth',0.5,'Color',0.5*[1,1,1]);
set(hs,'MarkerSize',20);
set(hs(2),'Color',niceblue);
set(hs(3),'Color',nicered);
set(hs(4),'Color',nicegreen);
ylabel('fraction correct');
hl = legend('MR','BP','BPSC','IsoRank',...
    'Location','SW');
legend('boxoff');
% legend('lp solution', 'rounded lp solution', 'correct matches lp', ...
%     'correct matches bp',...
%     'Location','NW');
xlabel('expected degree of noise in L (p \cdot n)');

greygrid;
set(gca,'Color','none'), set(gcf,'Color','none');
set(gcf, 'InvertHardCopy','off');
print(gcf,'-depsc2','-cmyk',[outputname '-ncorr.eps']);


%%

hs=plot(...
    pns, mean(fvals./fvals_ref), 'k--', ...     % upper mr
    pns, mean(rvals./fvals_ref), 'k.--', ...    % rounded mr
    pns, mean(rvals_bp./fvals_ref), 'k*-', ...  % rounded bp
    pns, mean(rvals_iso./fvals_ref), 'kd-.', ...% rounded iso
    pns, mean(ncorr)/400, 'k.--', ...
    pns, mean(ncorr_bp)/400, 'k*-', ...
    pns, mean(ncorr_iso)/400, 'kd-.');
line(xlim,[1 1],'LineWidth',0.5,'Color',0.5*[1,1,1]);
set(hs([2,5]),'MarkerSize',20);
set(hs(5:7),'Color',0.5*[1,1,1],'LineWidth',0.7);
ylabel('fraction of reference matching');
hl = legend('MR-upper','MR','BP','IsoRank',...
    'Location','SW');
legend('boxoff');
% legend('lp solution', 'rounded lp solution', 'correct matches lp', ...
%     'correct matches bp',...
%     'Location','NW');
xlabel('expected degree of noise in L (k^2 \cdot n)');

print(gcf,'-depsc2','-cmyk',[outputname '-all-bw.eps']);