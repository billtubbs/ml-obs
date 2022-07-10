function varargout = get_obs_vars_vecs(obs)
% varargout = get_obs_vars_vecs(obs)
% Returns data vectors containing values for all the
% time-varying variables of the observer. This is used
% in S-functions for storing the variable data in the
% Dwork memory objects allocated by Simulink.
%
% The number of data vectors returned depends on the
% type of observer. Some observers have integer
% variables which need to be stored in a separate
% Dwork vector.
%
% See function set_obs_vars_vecs for the reverse 
% operation - i.e. creating the value vectors from the
% observer struct.
%
% Examples
% >> vec_double = get_obs_vars_vecs(obs1);
% >> [vec_double, vec_int16] = get_obs_vars_vecs(obs2);
%

% Get observer variables
vars = get_obs_vars(obs);

switch obs.type

    case {"KFSS", "LB"}  % steady-state filters

        vars_double = {vars.xkp1_est, vars.ykp1_est};

        % Convert dynamic variables to vectors
        vdata = make_data_vectors(vars_double);
        varargout{1} = cell2mat(vdata.vecs);

    case "KF"  % standard Kalman filter

        vars_double = {vars.xkp1_est, vars.ykp1_est, vars.P};

        % Convert dynamic variables to vectors
        vdata = make_data_vectors(vars_double);
        varargout{1} = cell2mat(vdata.vecs);

    case {"MKF_SF95", "MKF_SF", "MKF_SP"}  % multi-model Kalman filters

        vars_double = {vars.xkp1_est, vars.ykp1_est, vars.p_seq_g_Yk, ...
            vars.gamma_k, vars.xkp1_est_f, vars.ykp1_est_f, vars.P_f};

        % Convert dynamic variables to vectors
        vdata = make_data_vectors(vars_double);
        vdata_int16 = make_data_vectors(struct2cell(vars.int16)', 'int16');
        varargout{1} = cell2mat(vdata.vecs);
        varargout{2} = cell2mat(vdata_int16.vecs);

    otherwise
        error('Value error: observer type not recognized')

end