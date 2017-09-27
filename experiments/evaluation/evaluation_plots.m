%% Plot 

%prob = 'lcsh2wiki-small'; upperb = 323; jitter=0.002;
%prob = 'lcsh2wiki'; upperb = 17608; jitter=0.02;
%prob = 'dmela-scere'; upperb = 381; jitter=0.002;
prob = 'musm-homo'; upperb = 1087; jitter=0.002;

load(sprintf('%s-isorank',prob));
load(sprintf('%s-bp',prob));
load(sprintf('%s-scbp',prob));
load(sprintf('%s-mr',prob));
load(sprintf('%s-mwm',prob));

%% Fix the figures
addpath('../tools/');

get_colors;

set(0, ...
    'defaultaxesfontsize',   14, ...
    'defaultaxeslinewidth',  1.0, ...
    'defaultlinelinewidth',   1.4, ...
    'defaultpatchlinewidth',  0.7)

% setup for poster
% set(0,'defaultaxesfontsize',18,'defaultaxeslinewidth',1.5,...
% 'defaultlinelinewidth',2,'defaultpatchlinewidth',.7) 

%%
% allhist = [];
% for i=1:length(bpresults)
%     allhist = [allhist; bpresults(i).hist];
% end
% for i=1:length(bpresults)
%     allhist = [allhist; bpresults(i).hist];
% end

%%
ms=0.5;
figure(1); clf; hold on; 
for i=1:length(bpresults), h1=plot(bpresults(i).hist(:,1),bpresults(i).hist(:,3),'b.','Color',niceblue,'MarkerSize',ms); end
for i=1:length(scbpresults), h1a=plot(scbpresults(i).hist(:,1),scbpresults(i).hist(:,3),'r.','Color',nicered,'MarkerSize',ms); end
for i=1:length(isoresults), h2=plot(isoresults(i).hist(:,1),isoresults(i).hist(:,3),'g.','Color',nicegreen,'MarkerSize',ms); end
for i=1:length(mrresults), h3=plot(mrresults(i).hist(:,1),mrresults(i).hist(:,3),'k.','MarkerSize',ms); end
xlabel('Weight'); ylabel('Overlap'); 
legend([h1,h1a,h2,h3],'BP','SCBP','IsoRank','MR','Location','NW');
legend boxoff;
box on; 
axis square;
%%
print(gcf,'-depsc2','-cmyk',[prob '-allhist.eps']);

%%
figure(2); clf; hold on;
maxoverall=[0 0 0];
markersize=3;
myrand=@() (rand-0.5);
legends = {};
for i=1:length(bpresults)
    h = bpresults(i).hist;
    [ignore,maxind] = max(h(:,3)); maxh = h(maxind,:);    
    maxoverall=max(maxoverall,maxh);
    h1=plot(maxh(1)*(1+jitter*myrand()),maxh(3)*(1+jitter*myrand()),...
        'o','MarkerSize',markersize,'Color',niceblue);
end
legends{end+1} = 'BP'; % must be outside the for loop

for i=1:length(scbpresults)
    h = scbpresults(i).hist;
    [ignore,maxind] = max(h(:,3)); maxh = h(maxind,:);    
    maxoverall=max(maxoverall,maxh);
    h1a=plot(maxh(1)*(1+jitter*myrand()),maxh(3)*(1+jitter*myrand()),...
        'o','MarkerSize',markersize,'Color',nicered);
end
legends{end+1} = 'SCBP'; % must be outside the for loop

for i=1:length(isoresults)
    h = isoresults(i).hist;
    [ignore,maxind] = max(h(:,3)); maxh = h(maxind,:);    
    maxoverall=max(maxoverall,maxh);
    h2=plot(maxh(1)*(1+jitter*myrand()),maxh(3)*(1+jitter*myrand()),...
        'o','MarkerSize',markersize,'Color',nicegreen);
end
legends{end+1} = 'IsoRank'; % must be outside the for loop

for i=1:length(mrresults)
    h = mrresults(i).hist;
    [ignore,maxind] = max(h(:,3)); maxh = h(maxind,:);    
    maxoverall=max(maxoverall,maxh);
    h3= plot(maxh(1)*(1+jitter*myrand()),maxh(3)*(1+jitter*myrand()),'ko','MarkerSize',markersize);
end
legends{end+1} = sprintf('MR'); % must be outside the for loop

maxweight = mwmresults(1).hist(1);
xlabel('Weight'); ylabel('Overlap'); 
% lcsh2wiki-small
% ylim
%xlim([50000,65000]); 
plot(maxweight,-1,'k.'); % include the max-weight area
yl = ylim;
xl = xlim;
ylim([0 yl(2)]);
xlim([0 xl(2)]);
yticks = get(gca,'YTick');
yticks(end+1) = maxoverall(3);
yticks = sort(yticks);
set(gca,'YTick',yticks);
ht=text(0.999*maxweight,0.3*yl(2),sprintf('max weight\n%g',maxweight));
set(ht,'HorizontalAlignment','right');
set(ht,'VerticalAlignment','top');
ht=text(0.5*xl(2), 0.99*upperb, sprintf('overlap upper bound\n%i',upperb));
set(ht,'HorizontalAlignment','right');
set(ht,'VerticalAlignment','top');
% ht=text(xl(2),upperb,sprintf('overlap upper bound\n%i',upperb));
% set(ht,'HorizontalAlignment','right');
% set(ht,'VerticalAlignment','top');

line(maxweight*[1,1],ylim,'Color',0.5*[1,1,1]);
line(xlim,upperb*[1,1],'Color',0.5*[1,1,1]);
box on;
legend([h1,h1a,h2,h3],legends{:},'Location','SW');
legend boxoff;
set(gca,'TickDir','out');
%axis square;
%%
% if necessary for decimals
set(gca, 'YTickLabel', num2str(transpose(get(gca, 'YTick')), '%1.0f') );
set(gca, 'XTickLabel', num2str(transpose(get(gca, 'XTick')), '%1.0f') );
%%
greygrid;
set(gca,'Color','none'), set(gcf,'Color','none');
set(gcf, 'InvertHardCopy','off');

print(gcf,sprintf('%s-summary.eps',prob),'-cmyk','-depsc2');
