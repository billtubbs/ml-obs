function T = objects2tablerow(vars)
% T = objects2tablerow(vars) converts a containers.Map of 
% variables (vars) into a one-row table T using the keys
% of the Map as column headings. Handles objects of the 
% following classes:
%  - double (scalars or arrays)
%  - char
%  - string
%  - string array
%  - struct
%  - cell arrays
%
% Example:
% >> sibs = {"Peter", "Fred"};
% >> vars = containers.Map({'Name', 'Age', 'Siblings'}, {"Amy", 5, sibs});
% >> T = objects2tablerow(vars)
% 
% T =
% 
%   1×4 table
% 
%     Age    Name     Siblings_1    Siblings_2
%     ___    _____    __________    __________
% 
%      5     "Amy"     "Peter"        "Fred"  
% 

    var_names = vars.keys;
    n_vars = numel(var_names);
    sections = {};
    n = 0;
    for i=1:n_vars
        switch class(vars(var_names{i}))
            
            case 'double'
                n = n + 1;
                if numel(vars(var_names{i})) == 1
                    sections{n} = table(vars(var_names{i}), 'VariableNames', var_names(i));
                else
                    sections{n} = array2tablerow(vars(var_names{i}), var_names{i});
                end
            
            case 'char'
                n = n + 1;
                sections{n} = table({vars(var_names{i})}, 'VariableNames', var_names(i));
            
            case 'string'
                n = n + 1;
                if numel(vars(var_names{i})) == 1
                    sections{n} = table(vars(var_names{i}), 'VariableNames', var_names(i));
                else
                    sections{n} = array2tablerow(vars(var_names{i}), var_names{i});
                end
            
            case 'struct'
                n = n + 1;
                keys = fieldnames(vars(var_names{i}));
                n_items = numel(keys);
                for j = 1:n_items
                    keys{j} = sprintf('%s_%s', var_names{i}, keys{j});
                end
                sections{n} = objects2tablerow(containers.Map(keys, struct2cell(vars(var_names{i}))));
            
            case 'cell'
                n = n + 1;
                n_items = numel(vars(var_names{i}));
                keys = cell(1, n_items);
                for j = 1:n_items
                    keys{j} = sprintf('%s_%d', var_names{i} , j);
                end
                sections{n} = objects2tablerow(containers.Map(keys, vars(var_names{i})));
        
        end
    end
    T = [sections{:}];
end