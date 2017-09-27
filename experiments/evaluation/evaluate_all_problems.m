%% Evaluate all datasets
% This script is just a wrapper around the evaluate_codes.m function to run
% it for all datasets automatically

problems = {'lcsh2wiki-small','musm-homo','dmela-scere','lcsh2wiki'};

for probi=1:length(problems)
    prob = problems{probi};
    evaluate_codes
end