% load full example
S = readSMAT('../../private_data/lcsubj2wikipedia-qp-squares.smat');
Ldata = load('../../private_data/lcsubj2wikipedia-qp-squares.edges');
li = Ldata(:,1)+1; 
lj = Ldata(:,2)+1;
w = Ldata(:,3);
