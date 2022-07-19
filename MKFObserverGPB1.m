% Multi-model Kalman Filter class definition
%
% obs = MKFObserverGPB1(A,B,C,Ts,P0,Q,R,T,label,x0,p_seq_g_Yk_init) 
% 
% Class for simulating the generalised pseudo-Bayes multi-
% model Kalman filter for state estimation of Markov jump
% linear systems. This is the filtering form of the 
% MKFObserver, which produces posterior estimates of the 
% states and outputs at the current time instant given the
% data at the current time:
%
%  x_hat(k|k) : estimate of states at time k
%  y_hat(k|k) : estimate of outputs at time k
%
% Arguments:
%   A, B, C : Cell arrays containing discrete-time system
%       matrices for each switching system modelled.
%   Ts : Sample period.
%   P0 : Initial covariance matrix of the state estimates
%       (same for each filter).
%   Q : Cell array of process noise covariance matrices for
%       each switching system.
%   R : Cell array of output measurement noise covariance
%       matrices for each switching system.
%   T : Transition probabity matrix of the Markov switching
%       process.
%   label : string name.
%   x0 : Initial state estimates (optional, default zeros).
%   p_seq_g_Yk_init : (optional, default uniform)
%       Initial prior hypothesis probabilities at time k-1.
%       If not specified, default is equal, i.e. uniform,
%       probability assigned to each hypothesis.
%


% TODO: Add these arguments:
% Prior estimates of the states and outputs at the next
% time instant given the data at the current time are
% also calculated:
%
%  xkp1_hat(k+1|k) : estimate of states at time k + 1
%  ykp1_hat(k+1|k) : estimate of outputs at time k + 1
%


classdef MKFObserverGPB1 < matlab.mixin.Copyable
    properties (SetAccess = immutable)
        Ts (1, 1) double {mustBeNonnegative}
        n (1, 1) double {mustBeInteger, mustBeNonnegative}
        nu (1, 1) double {mustBeInteger, mustBeNonnegative}
        ny (1, 1) double {mustBeInteger, mustBeNonnegative}
        nj (1, 1) double {mustBeInteger, mustBeNonnegative}
        f (1, 1) double {mustBeInteger, mustBeNonnegative}
        n_filt (1, 1) double {mustBeInteger, mustBeNonnegative}
    end
    properties
        A cell
        B cell
        C cell
        Q cell
        R cell
        seq cell
        T double
        label (1, 1) string
        P0 double
        x0 (:, 1) double
        p_seq_g_Yk_init double
        i (1, 1) {mustBeInteger, mustBeNonnegative}
        gamma_k double
        p_seq_g_Ykm1 double
        p_seq_g_Yk double
        p_yk_g_seq_Ykm1 double
        p_gammak_g_Ykm1 double
        p_gamma_k double
        filters  cell
        xk_est (:, 1) double
        Pk double
        yk_est (:, 1) double
        type (1, 1) string
    end
    methods
        function obj = MKFObserverGPB1(A,B,C,Ts,P0,Q,R,T,label,x0, ...
                p_seq_g_Yk_init)

            % System dimensions
            [n, nu, ny] = check_dimensions(A{1}, B{1}, C{1});

            % Number of switching systems
            nj = numel(A);

            % Number of filters required
            n_filt = nj;

            if nargin < 11
                % Initial values of prior conditional probabilities at 
                % k = -1. In absence of prior knowledge, assume all 
                % equally likely
                p_seq_g_Yk_init = ones(n_filt, 1) ./ double(n_filt);
            end
            if nargin < 10
                x0 = zeros(n,1);
            end
            obj.A = A;
            obj.B = B;
            obj.C = C;
            obj.Ts = Ts;
            obj.Q = Q;
            obj.R = R;
            obj.T = T;
            obj.label = label;
            obj.P0 = P0;
            obj.x0 = x0;

            % Model indicator values gamma(k) are static - there is
            % one filter for each model
            obj.gamma_k = (0:nj-1)';

            % Mode sequences (length 1 for GPB1)
            % For compatibility with other MKF observers.
            % TODO: Could eliminate seq, i and maybe f and make this a 
            % parent class of the other observers?
            obj.seq = num2cell(obj.gamma_k);

            % Sequence index starts at 0 but is 1 thereafter.
            obj.i = int16(0);

            % Fusion horizon length
            obj.f = 1;

            % Prior assumptions at initialization
            obj.p_seq_g_Yk_init = p_seq_g_Yk_init;

            % Check transition probability matrix
            % TODO: Is this necessary? If all possible hypotheses 
            % are not modelled (as is the case with some observers)
            % then perhaps T should not be whole?
            % assert(all(abs(sum(obj.T, 2) - 1) < 1e-15), "ValueError: T")
            assert(isequal(size(obj.T), [nj nj]), "ValueError: size(T)")

            % Check all other system matrix dimensions have same 
            % input/output dimensions and number of states.
            for j = 2:nj
                [n_j, nu_j, ny_j] = check_dimensions(A{j}, B{j}, C{j});
                assert(isequal([n_j, nu_j, ny_j], [n, nu, ny]), ...
                    "ValueError: size of A, B, C, and D")
            end

            % Store other useful variables
            obj.n = n;
            obj.nu = nu;
            obj.ny = ny;
            obj.nj = nj;
            obj.n_filt = n_filt;
            obj.type = "MKF_GPB1";

            % Initialize all variables
            obj.reset()

        end
        function reset(obj)
        % obj.reset()
        % Initialize all variables to initial values specified
        % when observer object was created.
        %

            % Sequence index does not change
            % TODO: Do we need it?
            obj.i = 0;

            % Initial values of prior conditional probabilities at k = -1 
            obj.p_seq_g_Yk = obj.p_seq_g_Yk_init;

            % Initialize state and output estimates
            % Note: At initialization at time k = 0, xkp1_est and
            % ykp1_est represent prior estimates of the states,
            % i.e. x_est(k|k-1) and y_est(k|k-1).
            obj.xk_est = obj.x0;
            obj.yk_est = obj.C{1} * obj.xk_est;

            % Initialize error covariance at k = 0
            obj.Pk = obj.P0;

            % At initialization at time k = 0, x_est(k|k)
            % and y_est(k|k) have not yet been computed.
            % TODO: Provide these
            %obj.xk_est = nan(obj.n, 1);
            %obj.yk_est = nan(obj.ny, 1);
            %obj.Pkp1 = obj.P0;

        end
        function update(obj, yk, uk)
        % obj.update(yk, uk)
        % updates the multi-model Kalman filter and calculates the
        % estimates of the states and output at the next sample
        % time.
        %
        % Arguments:
        %   obs : struct containing the multi-model Kalman filter
        %       variables (see function mkf_filter).
        %   uk : vector (nu, 1) of system inputs at the current 
        %       sample time.
        %   yk : vector (ny, 1) of system output measurements
        %       at the current sample time.
        %

            % Sequence index does not change
            % TODO: Do we need it?
            obj.i = 1;

            assert(isequal(size(uk), [obj.nu 1]), "ValueError: size(uk)")
            assert(isequal(size(yk), [obj.ny 1]), "ValueError: size(yk)")

            % Update all variables based on current measurement
            [obj.xk_est,obj.yk_est,obj.Pk,obj.p_seq_g_Yk] = ...
                GPB1_update(obj.A,obj.B,obj.C,obj.Q,obj.R,obj.T, ...
                obj.xk_est,obj.Pk,yk,uk,obj.p_seq_g_Yk);

        end

    end

end
