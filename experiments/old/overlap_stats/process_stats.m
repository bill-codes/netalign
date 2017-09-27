%% Load data
setup_experiments
results=[];
for ei=1:length(experiments)
    label = experiments(ei).label;
    results(end+1).label = experiments(ei).label;
    results(end).alpha = experiments(ei).alpha;
    results(end).beta = experiments(ei).beta;
    results(end).isorank = load(sprintf('isorank-%s.hist',label));
    results(end).bp_a = load(sprintf('bp-a-%s.hist',label));
    results(end).bp_b = load(sprintf('bp-b-%s.hist',label));
end
%%

%% Plot data
set(0,'defaultaxesfontsize',12,'defaultaxeslinewidth',.7,...
      'defaultlinelinewidth',.8,'defaultpatchlinewidth',.7) ;
figure(1); clf; hold all; 
colors = get(0,'DefaultAxesColorOrder');
curplot= 0; curlegend= {};
for ri=1:length(results)
    label=results(ri).label;
    if results(ri).alpha>results(ri).beta, continue; end
    curplot=curplot+1;
    hista= results(ri).bp_a;
    histb= results(ri).bp_b;
    histiso= results(ri).isorank;
    
    plot(histiso(:,4),histiso(:,5),'o','Color',colors(curplot,:)); 
    curlegend{end+1}=sprintf('isorank-%s',label);
    plot(hista(:,3),hista(:,4),'*','Color',colors(curplot,:)); 
    plot(histb(:,3),histb(:,4),'x','Color',colors(curplot,:)); 
    curlegend{end+1}=sprintf('bp-a-%s',label);
    curlegend{end+1}=sprintf('bp-b-%s',label);
    box on;
end
xlabel('Cardinality'); ylabel('Overlap');
lh=legend(curlegend{:},'Location','NW'); set(lh,'Box','off');
print(gcf,'qp-matching.eps','-depsc2');

figure(2); clf; hold all; 
colors = get(0,'DefaultAxesColorOrder');
curplot=0; curplot= 0; curlegend= {};
for ri=1:length(results)
    label=results(ri).label;
    if results(ri).alpha<results(ri).beta, continue; end
    curplot=curplot+1;
    hista= results(ri).bp_a;
    histb= results(ri).bp_b;
    histiso= results(ri).isorank;
    
    plot(histiso(:,4),histiso(:,5),'o','Color',colors(curplot,:)); 
    curlegend{end+1}=sprintf('isorank-%s',label);
    plot(hista(:,3),hista(:,4),'*','Color',colors(curplot,:)); 
    plot(histb(:,3),histb(:,4),'x','Color',colors(curplot,:)); 
    curlegend{end+1}=sprintf('bp-a-%s',label);
    curlegend{end+1}=sprintf('bp-b-%s',label);
    box on;
end
xlabel('Cardinality'); ylabel('Overlap');
lh=legend(curlegend{:},'Location','NW'); set(lh,'Box','off');
print(gcf,'qp-overlap.eps','-depsc2');
