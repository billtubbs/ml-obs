function vdata = make_data_vectors(varargin)
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