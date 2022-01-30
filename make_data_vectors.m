function vdata = make_data_vectors(varargin)
% vdata = make_data_vectors(varargin)
% Construct a set of vectors and associated meta data for
% the numeric variables passed to the function. This is useful
% for implementing S-functions in Simulink blocks which
% require all working data to be stored as one or more vectors.
%
% Example:
% >> vdata = make_data_vectors(1, [2 3; 4 5], {6, [7 8], 9})
% 
% vdata = 
% 
%   struct with fields:
% 
%      vecs: {[1]  [2 4 3 5]  [6 7 8 9]}
%     types: {'double'  'double'  {1×3 cell}}
%      dims: {[1 1]  [2 2]  {1×3 cell}}
% 
    vecs = cell(1, nargin);
    types = cell(1, nargin);
    dims = cell(1, nargin);
    for i = 1:nargin
        arg = varargin{i};
        switch class(arg)
            case 'double'
                vecs{i} = reshape(arg, 1, []);
                types{i} = 'double';
                dims{i} = size(arg);
            case 'cell'
                vd = make_data_vectors(arg{:});
                vecs{i} = cell2mat(vd.vecs);
                types{i} = reshape(vd.types, size(arg));
                dims{i} = reshape(vd.dims, size(arg));
            otherwise
                error("TypeError: invalid type.")
        end
    end
    vdata.vecs = vecs;
    vdata.types = types;
    vdata.dims = dims;
end