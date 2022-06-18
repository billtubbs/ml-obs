function obs = set_obs_vars_vecs(obs, varargin)
% obs = set_obs_vars_vecs(obs, varargin)
% Takes a certain number of vectors containing values
% for variables and sets those variables in the observer 
% struct. This is used in S-functions for 'unpacking'
% the variable data from the Dwork memory objects
% allocated by Simulink.
%
% See function get_obs_vars_vecs for the reverse 
% operation - i.e. creating the value vectors from the 
% observer struct.
%
% Example
% >> obs = set_obs_vars_vecs(obs, vec_double, vec_int16);
%
% In this example, vec_double contains the appropriate number
% of real double values and vec_int16 contains the
% appropriate number of integer values to specify all the
% variables in obs.
%

    n = obs.n;
    ny = obs.ny;

    switch obs.type

        case {'KFSS', 'LB'}  % steady-state filters

            assert(nargin == 2)
            vec_double = varargin{1};

            % double vector contents:
            % 1. xkp1_est : size(n, 1)
            % 2. ykp1_est : size(ny, 1)

            vdata.types = {'double', 'double'};
            vdata.dims = {[n 1], [ny 1]};
            vdata.n_els = {n, ny};

            % Add variables data
            vdata.vecs = mat2cell(vec_double', 1, cell2mat(vdata.n_els));

            % Unpack data vectors - doubles
            vars_double = unpack_data_vectors(vdata);
            vars = struct();
            vars.xkp1_est = vars_double{1};
            vars.ykp1_est = vars_double{2};

            % Set all dynamic variables
            obs = set_obs_vars(obs, vars);

        case 'KF'  % standard Kalman filter

            assert(nargin == 2)
            vec_double = varargin{1};

            % double vector contents:
            % 1. xkp1_est : size(n, 1)
            % 2. ykp1_est : size(ny, 1)
            % 3. P : size(n, n)

            vdata.types = {'double', 'double', 'double'};
            vdata.dims = {[n 1], [ny 1], [n n]};
            vdata.n_els = {n, ny, n*n};

            % Add variables data
            vdata.vecs = mat2cell(vec_double', 1, cell2mat(vdata.n_els));

            % Unpack data vectors - doubles
            vars_double = unpack_data_vectors(vdata);
            vars = struct();
            vars.xkp1_est = vars_double{1};
            vars.ykp1_est = vars_double{2};
            vars.P = vars_double{3};

            % Set all dynamic variables
            obs = set_obs_vars(obs, vars);

        case {'MKF', 'MKF_SF'}

            assert(nargin == 3)
            vec_double = varargin{1};
            vec_int16 = varargin{2};

            n_filt = obs.n_filt;
            f = obs.f;

            % double vector contents:
            % 1. xkp1_est : size(n, 1)
            % 2. ykp1_est : size(ny, 1)
            % 3. p_seq_g_Yk : size(n_filt, 1)
            % 4. xkp1_est for each KF : cell(1, n_filt)
            % 5. ykp1_est for each KF : cell(1, n_filt)
            % 6. P for each KF : cell(1, n_filt)

            % int16 vector contents:
            % 1. i : size(1, 2)
            % 2. i_next : size(1, 2)

            % Static data to unpack vectors
            % TODO: Could this be stored in the block somewhere
            %       other than Dwork memory.
            vdata.types = {'double', 'double', 'double', 'double', ...
                repmat({'double'}, 1, n_filt), ...
                repmat({'double'}, 1, n_filt), ...
                repmat({'double'}, 1, n_filt)};
            vdata.dims = {[n 1], [ny 1], [n_filt 1], [n_filt 1], ...
                repmat({[n 1]}, 1, n_filt), ...
                repmat({[ny 1]}, 1, n_filt), ...
                repmat({[n n]}, 1, n_filt)};
            vdata.n_els = {n, ny, n_filt, n_filt, n_filt*n, n_filt*ny, ...
                n_filt*n*n};
            vdata_int16.types = {'int16', 'int16'};
            vdata_int16.dims = {[1 2], [1 2]};
            vdata_int16.n_els = {2, 2};

            % Add variables data
            vdata.vecs = mat2cell(vec_double', 1, cell2mat(vdata.n_els));
            vdata_int16.vecs = mat2cell(vec_int16', 1, cell2mat(vdata_int16.n_els));

            % Unpack data vectors - doubles
            % TODO: Maybe this is unnecessary, just set the obs
            % attributes directly.
            vars_double = unpack_data_vectors(vdata);
            vars = struct();
            vars.xkp1_est = vars_double{1};
            vars.ykp1_est = vars_double{2};
            vars.p_seq_g_Yk = vars_double{3};
            vars.gamma_k = vars_double{4};
            vars.xkp1_est_f = vars_double{5};
            vars.ykp1_est_f = vars_double{6};
            vars.P_f = vars_double{7};

            % Unpack data vectors - integers
            vars_int16 = unpack_data_vectors(vdata_int16);
            vars.int16.i = vars_int16{1};
            vars.int16.i_next = vars_int16{2};

            % Set all dynamic variables
            obs = set_obs_vars(obs, vars);

        case 'MKF_SP'

            assert(nargin == 3)
            vec_double = varargin{1};
            vec_int16 = varargin{2};

            n_filt = obs.n_filt;
            f = obs.f;

            % double vector contents:
            % 1. xkp1_est : size(n, 1)
            % 2. ykp1_est : size(ny, 1)
            % 3. p_seq_g_Yk : size(n_filt, 1)
            % 4. gamma_k : size(n_filt, 1)
            % 5. xkp1_est for each KF : cell(1, n_filt)
            % 6. ykp1_est for each KF : cell(1, n_filt)
            % 7. P matrix for each filter : cell(1, n_filt)

            % int16 vector contents:
            % 1. i : size(1, 2)
            % 2. i_next : size(1, 2)
            % 3. f_main : size(1, n_main)
            % 4. f_hold : size(1, n_hold)
            % 5. seq for each KF : cell(n_filt, 1)

            % Static data to unpack vectors
            vdata.types = {'double', 'double', 'double', 'double', ...
                repmat({'double'}, 1, n_filt), ...
                repmat({'double'}, 1, n_filt), ...
                repmat({'double'}, 1, n_filt)};
            vdata.dims = {[n 1], [ny 1], [n_filt 1], [n_filt 1], ...
                repmat({[n 1]}, 1, n_filt), ...
                repmat({[ny 1]}, 1, n_filt), ...
                repmat({[n n]}, 1, n_filt)};
            vdata.n_els = {n, ny, n_filt, n_filt, n_filt*n, n_filt*ny, ...
                n_filt*n*n};
            vdata_int16.types = {'int16', 'int16', 'int16', 'int16', ...
                repmat({'int16'}, n_filt, 1)};
            vdata_int16.dims = {[1 2], [1 2], [1 obs.n_main], [1 obs.n_hold], ...
                repmat({[1 f]}, n_filt, 1)};
            vdata_int16.n_els = {int16(2), int16(2), obs.n_main, obs.n_hold, ...
                int16(n_filt)*f};

            % Add variables data
            vdata.vecs = mat2cell(vec_double', 1, cell2mat(vdata.n_els));
            vdata_int16.vecs = mat2cell(vec_int16', 1, cell2mat(vdata_int16.n_els));

            % Unpack data vectors - doubles
            vars_double = unpack_data_vectors(vdata);
            vars = struct();
            vars.xkp1_est = vars_double{1};
            vars.ykp1_est = vars_double{2};
            vars.p_seq_g_Yk = vars_double{3};
            vars.gamma_k = vars_double{4};
            vars.xkp1_est_f = vars_double{5};
            vars.ykp1_est_f = vars_double{6};
            vars.P_f = vars_double{7};

            % Unpack data vectors - integers
            vars_int16 = unpack_data_vectors(vdata_int16);
            vars.int16.i = vars_int16{1};
            vars.int16.i_next = vars_int16{2};
            vars.int16.f_main = vars_int16{3};
            vars.int16.f_hold = vars_int16{4};
            vars.int16.seq = vars_int16{5};

            % Set all dynamic variables
            obs = set_obs_vars(obs, vars);

        otherwise
            error('Value error: observer type not recognized')

    end