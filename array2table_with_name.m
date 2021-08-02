function T = array2table_with_name(X, label, sep)
% T = array2table_with_name(X, label, sep) converts the 
% matrix X to a table with custom column labels such as
% {'A1', 'A2', ... } if label is 'A' or {'A_1', 'A_2', 
% ... } if sep (optional) is '_'.
    if nargin == 2
        sep = '';
    end
    assert(numel(size(X)) == 2)
    n_cols = size(X, 2);
    if n_cols > 1
        col_names = compose(strcat(label, sep, '%d'), 1:n_cols);
    else
        col_names = {label};
    end
    T = array2table(X, 'VariableNames', col_names);
end