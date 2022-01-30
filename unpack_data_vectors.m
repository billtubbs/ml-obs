function varargout = unpack_data_vectors(vdata)
% varargout = unpack_data_vectors(vdata)
% Unpack a set of variables using the data in vdata.
% vdata is a struct created by the function 
% make_data_vectors for storing numerical data as
% a set of vectors, which is useful for implementing
% S-functions in Simulink blocks.
%
% Example:
% >> vdata = make_data_vectors(1, [2 3; 4 5], {6, [7 8], 9});
% >> [a, b, c] = unpack_data_vectors(vdata)
% 
% a =
% 
%      1
% 
% 
% b =
% 
%      2     3
%      4     5
% 
% 
% c =
% 
%   1×3 cell array
% 
%     {[6]}    {1×2 double}    {[9]}
% 
    varargout = cell(1, nargout);
    for i = 1:nargout
        vec = vdata.vecs{i};
        type = vdata.types{i};
        dim = vdata.dims{i};
        if iscell(type)
            vd.vecs = mat2cell(vec, 1, cellfun(@numel_recursive, dim));
            vd.types = type;
            vd.dims = dim;
            args = cell(1, numel(vd.vecs));
            [args{:}] = unpack_data_vectors(vd);
            varargout{i} = args;
        elseif isequal(type, 'double')
            varargout{i} = reshape(vec, dim);
        else
            error("TypeError: invalid type.")
        end
    end
end