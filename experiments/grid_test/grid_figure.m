%% Generate a figure demonstrating the grid setup

%%
addpath('../misc/gaimc'); % get bfs function
addpath('../../matlab'); % get bfs function

%%
k=8;
q=0.25;
d=1;
p=0.5/k^2;
[A,B,L,xy,G,Lstart,Lrand] = align_grid_test_data(k,q,p,d);
G = triu(spones(G)); A = triu(spones(A)); B = triu(spones(B)); 
L = spones(L); 
Lstart = spones(Lstart); 
Lrand = spones(Lrand);
Ldist = L - Lstart - Lrand;

%%
clf; hold on;
z1 = 0; z2=0.5;
[Gx,Gy]=gplot(G,xy); Gz1=z1*ones(length(Gx),1); Gz2=z2*ones(length(Gx),1);

plot3(Gx,Gy,Gz1,'k.-','MarkerSize',18,'LineWidth',1.2); 
plot3(Gx,Gy,Gz2,'k.-','MarkerSize',18,'LineWidth',1.2);
camorbit(-20,65);

% draw the extra edges
[Gx,Gy]=gplot(A-G,xy); Gz=z1*ones(length(Gx),1); plot3(Gx,Gy,Gz,'r-');
[Gx,Gy]=gplot(B-G,xy); Gz=z2*ones(length(Gx),1); plot3(Gx,Gy,Gz,'r-');
axis off;

[Gx,Gy]=gplot(Lstart,xy); 
Gz=repmat([z1;z2;NaN],1,length(Gx)/3); Gz=Gz(:);
plot3(Gx,Gy,Gz,'k-','LineWidth',0.5,'Color',0.5*[1,1,1]);

set(gcf,'InvertHardCopy','off');
set(gca,'Color','none'); set(gcf,'Color','none');
print(gcf,'grid-fig-match.eps','-depsc2');

%%
clf; hold on;
z1 = 0; z2=0.5;
[Gx,Gy]=gplot(G,xy); Gz1=z1*ones(length(Gx),1); Gz2=z2*ones(length(Gx),1);

plot3(Gx,Gy,Gz1,'k.-','MarkerSize',18,'LineWidth',1.2); 
plot3(Gx,Gy,Gz2,'k.-','MarkerSize',18,'LineWidth',1.2);
camorbit(-20,65);

% draw the extra edges
[Gx,Gy]=gplot(A-G,xy); Gz=z1*ones(length(Gx),1); plot3(Gx,Gy,Gz,'r-');
[Gx,Gy]=gplot(B-G,xy); Gz=z2*ones(length(Gx),1); plot3(Gx,Gy,Gz,'r-');
axis off;

[Gx,Gy]=gplot(Lrand,xy); 
Gz=repmat([z1;z2;NaN],1,length(Gx)/3); Gz=Gz(:);
plot3(Gx,Gy,Gz,'k-','LineWidth',0.5,'Color',0.5*[0,0,1]);

set(gcf,'InvertHardCopy','off');
set(gca,'Color','none'); set(gcf,'Color','none');
print(gcf,'grid-fig-rand.eps','-depsc2');
%%
clf; hold on;
z1 = 0; z2=0.5;
[Gx,Gy]=gplot(G,xy); Gz1=z1*ones(length(Gx),1); Gz2=z2*ones(length(Gx),1);

plot3(Gx,Gy,Gz1,'k.-','MarkerSize',18,'LineWidth',1.2); 
plot3(Gx,Gy,Gz2,'k.-','MarkerSize',18,'LineWidth',1.2);
camorbit(-20,65);

% draw the extra edges
[Gx,Gy]=gplot(A-G,xy); Gz=z1*ones(length(Gx),1); plot3(Gx,Gy,Gz,'r-');
[Gx,Gy]=gplot(B-G,xy); Gz=z2*ones(length(Gx),1); plot3(Gx,Gy,Gz,'r-');
axis off;


[Gx,Gy]=gplot(Ldist,xy); 
Gz=repmat([z1;z2;NaN],1,length(Gx)/3); Gz=Gz(:);
plot3(Gx,Gy,Gz,'k-','LineWidth',0.5,'Color',0.5*[1,0,0]);

set(gcf,'InvertHardCopy','off');
set(gca,'Color','none'); set(gcf,'Color','none');
print(gcf,'grid-fig-dist.eps','-depsc2');

