function obs = set_obs_vars_vecs(obs, varargin)

switch obs.type

    case 'MKF_AFMM'

        assert(nargin == 3)
        vec_double = varargin{1};
        vec_int16 = varargin{2};

        % Static data to unpack vectors
        vdata.types = {'double', 'double', 'double', ...
            {'double', 'double', 'double', 'double', 'double'}, ...
            {'double', 'double', 'double', 'double', 'double'}, ...
            {'double', 'double', 'double', 'double', 'double'}};
        vdata.dims = {[2 1], [1 1], [5 1], {[2 1], [2 1], [2 1], [2 1], [2 1]}, ...
            {[1 1], [1 1], [1 1], [1 1], [1 1]}, ...
            {[2 2], [2 2], [2 2], [2 2], [2 2]}};
        vdata.n_els = {2, 1, 5, 10, 5, 20};
        vdata_int16.types = {'int16', 'int16', 'int16', 'int16', 'int16', {'int16'; 'int16'; 'int16'; 'int16'; 'int16'}};
        vdata_int16.dims = {[1 2], [1 2], [1 3], [1 2], [1 4], {[1 10]; [1 10]; [1 10]; [1 10]; [1 10]}};
        vdata_int16.n_els = {2, 2, 3, 2, 4, 50};

        % Add variables data
        vdata.vecs = mat2cell(vec_double', 1, cell2mat(vdata.n_els));
        vdata_int16.vecs = mat2cell(vec_int16', 1, cell2mat(vdata_int16.n_els));

        % Unpack data vectors - doubles
        vars_double = unpack_data_vectors(vdata);
        vars = struct();
        vars.xkp1_est = vars_double{1};
        vars.ykp1_est = vars_double{2};
        vars.p_seq_g_Yk = vars_double{3};
        vars.xkp1_est_f = vars_double{4};
        vars.ykp1_est_f = vars_double{5};
        vars.P_f = vars_double{6};

        % Unpack data vectors - integers
        vars_int16 = unpack_data_vectors(vdata_int16, 'int16');
        vars.int16.i = vars_int16{1};
        vars.int16.i_next = vars_int16{2};
        vars.int16.f_main = vars_int16{3};
        vars.int16.f_hold = vars_int16{4};
        vars.int16.f_unused = vars_int16{5};
        vars.int16.seq = vars_int16{6};

        % Set all dynamic variables
        obs = set_obs_vars(obs, vars);

    otherwise
        error('Value error: observer type not recognized')

end