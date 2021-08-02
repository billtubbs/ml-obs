function T = array2tablerow(X, label)
    if nargin == 1
        label = inputname(1);
    end
    if isscalar(X)
        T = array2table(X, 'VariableNames', {label});
    else
        switch numel(size(X))

            case 2
                [R, C] = ind2sub(size(X), 1:numel(X));
                ind = [R' C'];
                n = size(ind, 1);
                col_names = cell(1, n);
                if ismember(1, size(X))
                    j = find(size(X) ~= 1);
                    for i=1:n
                        col_names{i} = sprintf('%s_%d', label, ind(i, j));
                    end
                else
                    for i=1:n
                        col_names{i} = sprintf('%s_%d_%d', label, ind(i,:));
                    end
                end

            case 3
                [R, C, D] = ind2sub(size(X), 1:numel(X));
                ind = [R' C' D'];
                n = size(ind, 1);
                col_names = cell(1, n);
                for i=1:n
                    col_names{i} = sprintf('%s_%d_%d_%d', label, ind(i,:));
                end

        end
        T = array2table(reshape(X,1,[]), 'VariableNames', col_names);
        
    end
end