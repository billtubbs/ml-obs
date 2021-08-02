function S = combinations_lte(n, m)
% S = combinations_lt(n, m) returns a cell array of
% vectors of width n representing all combinations of
% sequences with no more than m values set to 1 and
% all other values 0.
%
% Example
%
% >> S = combinations_lte(3, 2);
% >> cell2mat(S)
% 
% ans =
% 
%      0     0     0
%      1     0     0
%      0     1     0
%      0     0     1
%      1     1     0
%      1     0     1
%      0     1     1
% 

    C = cell(m+1, 1);
    for i = 0:m
        C{i+1} = combinations(n, i);
    end
    S = vertcat(C{:});
end