%% Output a table from the evaluation experiments
% The table shows:
%
% Method | Problem | Best overlap | Parameters


methods = {'mwm','iso','bp','scbp','mr'};
probs = {'musm-homo','dmela-scere','lcsh2wiki-small','lcsh2wiki'};
upperb = [1087,381,323,17608];

mresults = cell(size(methods));
for pi=1:length(probs)
    prob = probs{pi};
    load(sprintf('%s-isorank',prob));
    load(sprintf('%s-bp',prob));
    load(sprintf('%s-scbp',prob));
    load(sprintf('%s-mr',prob));
    load(sprintf('%s-mwm',prob));
    for mi=1:length(methods)
        method = methods{mi};
        cresults = eval(sprintf('%sresults',method));
        maxind = 1; [maxover maxoveriter]=max(cresults(1).hist(:,3));
        for i=2:length(cresults)
            if maxover < max(cresults(i).hist(:,3))
                maxind = i; [maxover maxoveriter]=max(cresults(i).hist(:,3));
            end
        end
        mresults{mi}{pi} = {...
            maxover, maxoveriter, cresults(maxind).str,  cresults(maxind)};
    end
end

for mi=1:length(methods)
    mstr = methods{mi};
    for pi=1:length(probs);
        pstr = probs{pi};
        s = mresults{mi}{pi}; % summary
        sr = s{4};
        fprintf('%5s ',mstr); mstr = ''; % only show mstr once!
        fprintf('& %15s ',pstr); pstr = ''; % only show pstr once!
        fprintf('& %7i ',s{1}); % show max overlap
        fprintf('& %5.1f\\%% ', 100*s{1}/upperb(pi)); % show overlap percent of upperbound
        if isfield(sr,'time')
            fprintf('& %7.1f ',(s{2}/length(sr.hist))*sr.time);
            fprintf('& %7.1f ',sr.time);
        else
            fprintf('& %7s ','---');
            fprintf('& %7s ','---');
        end
        %fprintf('& %5i ', s{2}); % iteration number
        switch methods{mi}
        case 'mwm'
            fprintf('& '); % nothing extra
        case 'iso'
            fprintf('& gamma %4.2f; rounding %i ', sr.alpha, sr.rtype);
        case 'bp'
            fprintf('& alpha %2i; beta %2i; gamma %5.3f; dtype %i ',...
                sr.obj(1), sr.obj(2), sr.gamma, sr.dtype);
        case 'scbp'
            fprintf('& alpha %2i; beta %2i; gamma %5.3f; dtype %i ',...
                sr.obj(1), sr.obj(2), sr.gamma, sr.dtype);
        case 'mr'
            fprintf('& alpha %2i; beta %2i; gamma %5.3f; step %i ',...
                sr.obj(1), sr.obj(2), sr.gamma, sr.stepm);
        end
        fprintf('\\\\\n');
    end
end
        
