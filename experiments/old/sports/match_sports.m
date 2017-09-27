addpath('~/devextern/Clp-1.9.0/')
addpath('../../matlab');
%% 
% Data from Kathryn Pendings (2009-08-11)
SoConDavg = [
         0    6.0000         0         0         0   19.0000    3.5000   12.0000         0    1.0000    6.0000    5.0000
    4.0000         0    6.0000   11.0000         0   11.5000   27.0000   13.0000    8.6667   13.5000    7.0000    7.0000
    3.5000         0         0   11.5000   18.0000    2.0000   13.0000   15.0000   25.0000    6.5000   11.5000         0
    6.0000    9.0000         0         0    4.5000    8.0000    9.0000   13.0000    7.0000   14.6667   15.0000    2.0000
   17.0000   13.5000   15.0000    4.0000         0   13.5000   27.5000   32.0000   12.0000   21.0000   24.0000   20.0000
         0    1.0000    2.0000    6.0000         0         0    8.0000    8.5000         0         0    6.0000    7.0000
         0    2.0000         0         0         0         0         0   10.0000    6.0000    5.0000         0         0
         0    4.0000         0         0         0         0   16.0000         0    7.0000    9.0000         0   21.0000
   10.5000         0    9.0000         0         0   13.0000    8.0000   23.0000         0   14.0000   12.0000         0
    2.0000         0         0         0         0   12.5000         0    8.0000         0         0         0         0
    5.0000    9.0000         0    2.0000         0    7.0000   10.0000   14.0000    7.0000    9.0000         0   12.0000
    6.0000    3.0000    5.0000         0         0    3.0000   14.0000    7.0000   11.0000   12.0000   15.0000         0
    ];

teams = {
    'App State'
	'UT Chatt'
    'Citadel'
    'CofC'
    'Davidson'
    'Elon'
    'Furman'
    'GA Southern'
    'Samford'
    'UNC-G'
    'W. Carolina'
    'Wofford'};

scores = sum(SoConDavg,2);
%%
% Setup 
B = triu(ones(12),1)';
A = SoConDavg>0;
[S,w,li,lj] = netalign_setup_dir(A,B,ones(12,12));
[f,C,d,Si,Sj] = netalign_lp_prob_dir(S,w,0,1,li,lj,'tight');
[x,z,status] = clp([],-f,C,d,[],[],zeros(length(f),1),ones(length(f),1));
status
%%
xs = x(1:length(w));
[ma mb mi weight overlap] = mwmround(xs,S,w,li,lj); overlap
[ma mb mi weight overlap] = mwmround(S*xs,S,w,li,lj); overlap
Y = sparse(Si,Sj,x(length(w)+1:end),size(S,1),size(S,2));
%mwmround(Y*xs,S,w,li,lj)
[ma mb mi weight overlap] = mwmround(Y*xs,S,w,li,lj); overlap
[ma mb mi weight overlap] = mwmround(sum(Y,2),S,w,li,lj); overlap

%%
% look at ranks
[num2cell(mb) teams(ma)]
imb = 0; imb(mb) = 1:length(mb);
% look at scores
sum(SoConDavg(imb,imb),2)
%SoConDavg(mb(ma),mb(ma))

%% 
% print data for presentation
p2= fliplr(imb);
P = SoConDavg(p2,p2)

fprintf('    ');
for i=1:length(p2)
    ti = char(teams(p2(i)));
    fprintf('& %3s ',ti(1:3));
end
fprintf('\n');
for i=1:length(p2)
    ti = char(teams(p2(i)));
    fprintf('%3s ',ti(1:3));
    for j=1:length(p2)
        v = round(SoConDavg(p2(i),p2(j)));
        if v==0,
            fprintf('&  ');
        else
            fprintf('& %i ', v);
        end
    end
    fprintf('\\\\ \n');
end
        

%% 
% Look at some variations
temp = find(mb==1)
mb(temp) = mb(5);
mb(5) = 1;

%% Lemme look at some silly examples
n = 12;
p = randperm(n);
A = triu(ones(n),1);
B = triu(ones(n),1);
B(2,1) = 1;
B = B(p,p);
[S,w,li,lj] = netalign_setup_dir(A,B,ones(size(A,1),size(B,1)));
[f,C,d,Si,Sj] = netalign_lp_prob_dir(S,w,0,1,li,lj,'tight');
[x,z,status] = clp([],-f,C,d,[],[],zeros(length(f),1),ones(length(f),1));
status
xs = x(1:length(w));
[ma mb mi weight overlap] = mwmround(S*xs,S,w,li,lj); overlap
f'*x
