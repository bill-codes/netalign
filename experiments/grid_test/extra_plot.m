% an extra plot function, I'm not entirely sure what this one was for, 
% but I'm keeping it just in case it was used for the paper.  We'll do
% some cleanup for the journal repo

set(0,'defaultaxesfontsiz',12,'defaultaxeslinewidth',.8,...
'defaultlinelinewidth',.9,'defaultpatchlinewidth',.7) 

clf, %subplot('position',[.1 .4 .8 .5]); 
%load 'results-k-20'
% comparison: lp, llp, (upper); bp, iso, lp, llp (rounded, lower)
%             --  :                :, -., -- -,
hs=plot(...
    pns, max(fvals./fvals_ref), 'k--', ...     % upper mr
    pns, max(rvals./fvals_ref), 'k.--', ...    % rounded mr
    pns, max(rvals_bp./fvals_ref), 'k*-', ...  % rounded bp
    pns, max(rvals_iso./fvals_ref), 'kd-.', ...% rounded iso
    pns, max(ncorr)/400, 'k.--', ...
    pns, max(ncorr_bp)/400, 'k*-', ...
    pns, max(ncorr_iso)/400, 'kd-.');
line(xlim,[1 1],'LineWidth',0.5,'Color',0.5*[1,1,1]);
set(hs([2,5]),'MarkerSize',13);
set(hs(5:7),'Color',0.5*[1,1,1],'LineWidth',0.7);
ylabel('fraction of reference matching');
hl = legend('MR-upper','MR','BP','IsoRank',...
    'Location','SW');
legend('boxoff');
% legend('lp solution', 'rounded lp solution', 'correct matches lp', ...
%     'correct matches bp',...
%     'Location','NW');
xlabel('expected degree of noise in L (k^2 \cdot n)');
print(gcf,'-depsc2','-cmyk','grid-nollp.eps');
