set(0,'defaultaxesfontsiz',12,'defaultaxeslinewidth',.8,...
'defaultlinelinewidth',.9,'defaultpatchlinewidth',.7) 

clf, %subplot('position',[.1 .4 .8 .5]); 
%load 'results-k-20'
% comparison: lp, llp, (upper); bp, iso, lp, llp (rounded, lower)
%             --  :                :, -., -- -,
hs=plot(...
    pns, prctile(fvals./fvals_ref,50), 'k--', ...     % upper mr
    pns, prctile(rvals./fvals_ref,50), 'k.--', ...    % rounded mr
    pns, prctile(rvals_bp./fvals_ref,50), 'k*-', ...  % rounded bp
    pns, prctile(rvals_iso./fvals_ref,50), 'kd-.', ...% rounded iso
    pns, prctile(ncorr,50)/400, 'k.--', ...
    pns, prctile(ncorr_bp,50)/400, 'k*-', ...
    pns, prctile(ncorr_iso,50)/400, 'kd-.');
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