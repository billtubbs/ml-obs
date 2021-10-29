function obs = MEKF_observer(n,state_fcn,meas_fcn,params,u_meas,y_meas,dfdx,dhdx,Ts, ...
    P0,Q,R,seq,T,d,label,x0,y0)
% obs = MEKF_observer(n,state_fcn,meas_fcn,params,u_meas,y_meas,dfdx,dhdx,Ts, ...
%     P0,Q,R,seq,T,d,label,x0,y0)
% Creates a struct for simulating a multi-model extended 
% Kalman filter for state estimation of a non-linear 
% Markov jump system.
%
% Arguments:
%   n : Number of model states.
%   state_fcn, meas_fcn : cell arrays containing state 
%       transition functions and measurement functions for 
%       each switching system modelled.
%   params : cell array of cell arrays containing any 
%       additional parameters for each switching system that
%       should be passed to functions f, h, dfdx, and dhdx.
%   u_meas : array indicating which inputs are measured.
%   y_meas : array indicating which outputs are measured.
%   dfdx : cell array containing functions to generate the
%       Jacobian matrix of the state transition function for
%       each switching system.
%   dhdx : cell array containing functions to generate the
%       Jacobian matrix of the measurement function for each
%       switching system.
%   Ts : sample period.
%   P0 : cell array of initial covariance matrices of the 
%       state estimates for each filter.
%   Q : cell array of process noise covariance matrices for
%       each switching system.
%   R : cell array of output measurement noise covariance
%       matrices for each switching system.
%   seq : model indicator sequences for each filter (in rows).
%   T : transition probabity matrix of the Markov switching
%       process.
%   d : detection interval length in number of sample periods.
%   label : string name.
%   x0 : intial state estimates (optional).
%   y0 : intial output estimates (optional).
%

    % Number of switching systems
    nj = numel(state_fcn);
    assert(numel(meas_fcn) == nj)
    assert(numel(dfdx) == nj)
    assert(numel(dhdx) == nj)

    obs.n = n;
    if nargin == 16
        % Default initial state estimate
        x0 = zeros(n, 1);
    end
    ny = size(y_meas,1);
    if nargin < 18
        % Default initial state estimate
        y0 = zeros(ny, 1);
    end

    obs.state_fcn = state_fcn;
    obs.meas_fcn = meas_fcn;
    obs.params = params;
    obs.u_meas = u_meas;  % TODO implement these
    obs.y_meas = y_meas;  % 
    obs.dfdx = dfdx;
    obs.dhdx = dhdx;
    obs.Ts = Ts;
    obs.P0 = P0;
    obs.Q = Q;
    obs.R = R;
    obs.seq = seq;
    obs.T = T;
    obs.d = d;
    obs.label = label;

    % Check transition probability matrix
    assert(all(abs(sum(T, 2) - 1) < 1e-15), "ValueError: T")

    % Check matrices dimensions.
    for j = 1:nj
        assert(isequal(size(R{j}), [ny ny]))
        assert(isequal(size(Q{j}), [n n]))
    end

    % Number of filters required
    n_filt = size(seq, 1);

    % Fusion horizon length
    f = size(cell2mat(seq), 2);

    % Sequence index and counter for prob. updates
    % obs.i(1) is the sequence index (1 <= i(1) <= obs.f)
    % obs.i(2) is the counter for prob. updates (1 <= i(2) <= obs.d)
    obs.i = nan(1, 2);
    obs.i_next = int16([1 1]);

    % Initialize conditional probabilities: all equally likely
    obs.p_seq_g_Yk = ones(n_filt, 1) ./ n_filt;

    % Empty vectors to store values for filter calculations
    % p(y(k)|Gamma(k),Y(k-1))
    obs.p_yk_g_seq_Ykm1 = zeros(n_filt, 1);
    % Pr(gamma(k)|Y(k-1))
    obs.p_gammak_g_Ykm1 = zeros(n_filt, 1);
    % Pr(gamma(k))
    obs.p_gamma_k = zeros(n_filt, 1);
    % Pr(Gamma(k)|Y(k-1))
    obs.p_seq_g_Ykm1 = zeros(n_filt, 1);

    % Create multi-model observer
    obs.filters = cell(n_filt, 1);
    fmt = strcat('%s%0', char(string(strlength(string(n_filt)))), 'd');
    for i = 1:n_filt
        label_i = sprintf(fmt,label,i);
        % Initialize each filter with system #1
        obs.filters{i} = EKF_observer(n,state_fcn{1},meas_fcn{1}, ...
            params{1},u_meas,y_meas,dfdx{1},dhdx{1},Ts,P0{i},Q{1}, ...
            R{1},label_i,x0,y0);
    end

    % Initialize estimates
    obs.xkp1_est = x0;
    obs.ykp1_est = y0;
    obs.Pkp1 = P0{1};  % TODO: is this correct?

    % Save useful variables in struct
    obs.nj = nj;
    obs.nu = size(u_meas,1);
    obs.ny = ny;
    obs.n_filt = n_filt;
    obs.f = f;

end