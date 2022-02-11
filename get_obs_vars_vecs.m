function [vec_double, vec_int16] = get_obs_vars_vecs(obs)

% Get observer variables
vars = get_obs_vars(obs);

switch obs.type

    case {'MKF', 'MKF_RODD', 'MKF_AFMM'}

        vars_double = {vars.xkp1_est, vars.ykp1_est, vars.p_seq_g_Yk, ...
            vars.xkp1_est_f, vars.ykp1_est_f, vars.P_f};

        % Convert dynamic variables to vectors
        vdata = make_data_vectors(vars_double);
        vdata_int16 = make_data_vectors(struct2cell(vars.int16)', 'int16');
        vec_double = cell2mat(vdata.vecs);
        vec_int16 = cell2mat(vdata_int16.vecs);

    otherwise
        error('Value error: observer type not recognized')

end