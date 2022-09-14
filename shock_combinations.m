function S = shock_combinations(n, m)
% S = shock_combinations(n, m) returns a cell array of
% integer vectors of width n representing all
% combinations of sequences with m values set to
% 1 and all other values 0.
%
% Example:
% >> S = shock_combinations(3, 2);
% >> cell2mat(S)
% 
% ans =
% 
%   3×3 int16 matrix
% 
%    2   2   1
%    2   1   2
%    1   2   2
%    

    C = int16(nchoosek(1:n, m));
    S = repmat({int16(ones(1, n))}, size(C, 1), 1);
    for i = 1:size(C, 1)
        S{i}(1, C(i,:)) = 2;
    end

end