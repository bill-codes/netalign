%% Generate a figure demonstrating the grid setup

%%
addpath('../misc/gaimc'); % get bfs function
addpath('../../matlab'); % get bfs function
addpath('~/dev/matlab-bgl-4.0'); % get a faster bfs and a graph layout alg

%% compile the generation code
!g++ -O2 RandomPowerLaw.cpp

%% Set the parameters for the experiment

nrep = 1;
theta = 1.8;
q = 0.02;
d = 1;
n = 50;
p = 0.5;

cmd = sprintf('./a.out %i %f %f %f %i',n,theta,q,p/n,nrep);
system(cmd);
G = readSMAT(sprintf('A%i-orig.smat',0));
A = readSMAT(sprintf('A%i.smat',0));
B = readSMAT(sprintf('B%i.smat',0));
Ldata = load(sprintf('L%i.data',0));
!rm A*.smat 
!rm B*.smat 
!rm S*.smat 
!rm L*.data
Ldata(:,1)=Ldata(:,1)+1; Ldata(:,2)=Ldata(:,2)+1;
Ldata(end+1,:)=[n,n,0];
L = spconvert(Ldata); L = spones(L);
A = A|A';
B = B|B';
G = G|G';
Lstart = speye(n);
Lrand = L-Lstart;
%% Compute a layout
clf;
xy = kamada_kawai_spring_layout(double(A));
gplot(G,xy,'.-');

%% Load the real data

% load('powerlaw-fig-data');

%%
% plot just the graph
clf; hold on;
[Gx,Gy]=gplot(G,xy);
plot(Gx,Gy,'k.-','MarkerSize',18,'LineWidth',1.2);
axis off;

set(gca,'Color','none'); set(gcf,'Color','none');
set(gcf,'InvertHardCopy','off');
print(gcf,'powerlaw-fig-graph.eps','-depsc2');

%%
clf; hold on;
z1 = 0; z2=0.1;
[Gx,Gy]=gplot(G,xy); Gz1=z1*ones(length(Gx),1); Gz2=z2*ones(length(Gx),1);

plot3(Gx,Gy,Gz1,'k.-','MarkerSize',18,'LineWidth',1.2); 
plot3(Gx,Gy,Gz2,'k.-','MarkerSize',18,'LineWidth',1.2);
camorbit(-15,55);

% draw the extra edges
[Gx,Gy]=gplot(A-G,xy); Gz=z1*ones(length(Gx),1); plot3(Gx,Gy,Gz,'r-');
[Gx,Gy]=gplot(B-G,xy); Gz=z2*ones(length(Gx),1); plot3(Gx,Gy,Gz,'r-');
axis off;

[Gx,Gy]=gplot(Lstart,xy); 
Gz=repmat([z1;z2;NaN],1,length(Gx)/3); Gz=Gz(:);
plot3(Gx,Gy,Gz,'k-','LineWidth',0.5,'Color',0.5*[1,1,1]);

set(gca,'Color','none'); set(gcf,'Color','none');
set(gcf,'InvertHardCopy','off');
print(gcf,'powerlaw-fig-match.eps','-depsc2');

%%
clf; hold on;
z1 = 0; z2=0.5;
[Gx,Gy]=gplot(G,xy); Gz1=z1*ones(length(Gx),1); Gz2=z2*ones(length(Gx),1);

plot3(Gx,Gy,Gz1,'k.-','MarkerSize',18,'LineWidth',1.2); 
plot3(Gx,Gy,Gz2,'k.-','MarkerSize',18,'LineWidth',1.2);
camorbit(-15,55);

% draw the extra edges
[Gx,Gy]=gplot(A-G,xy); Gz=z1*ones(length(Gx),1); plot3(Gx,Gy,Gz,'r-');
[Gx,Gy]=gplot(B-G,xy); Gz=z2*ones(length(Gx),1); plot3(Gx,Gy,Gz,'r-');
axis off;

[Gx,Gy]=gplot(Lrand,xy); 
Gz=repmat([z1;z2;NaN],1,length(Gx)/3); Gz=Gz(:);
plot3(Gx,Gy,Gz,'k-','LineWidth',0.5,'Color',0.5*[0,0,1]);

set(gca,'Color','none'); set(gcf,'Color','none');
set(gcf,'InvertHardCopy','off');
print(gcf,'powerlaw-fig-rand.eps','-depsc2');


%%
% save the final data
save 'powerlaw-fig-data.mat' A B L G xy Lstart Lrand;