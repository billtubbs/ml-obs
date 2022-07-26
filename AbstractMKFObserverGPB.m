% Multi-model Kalman Filter class definition
%
% Abstract base class for the generalised pseudo-Bayes multi-
% model Kalman filters for state estimation of Markov 
% jump linear systems. The defining characteristic of
% these methods is that the number of independent
% Kalman filters matches the number of system models.
% I.e. obj.n_filt = obj.nj
%
% obs = AbstractMKFObserverGPB(A,B,C,Ts,P0,Q,R,T,label,x0, ...
%     p_seq_g_Yk_init)
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
%   n_filt : Number of independent filters to model.
%   label : string name.
%   x0 : Initial state estimates (optional, default zeros).
%   p_seq_g_Yk_init : (optional, default uniform)
%       Initial prior hypothesis probabilities at time k-1.
%       If not specified, default is equal, i.e. uniform,
%       probability assigned to each hypothesis.
%

classdef (Abstract) AbstractMKFObserverGPB < matlab.mixin.Copyable
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
        T double
        label (1, 1) string
        P0 double
        x0 (:, 1) double
        gamma_k double
        p_seq_g_Yk_init double
        p_seq_g_Ykm1 double
        p_seq_g_Yk double
        p_yk_g_seq_Ykm1 double
        p_gammak_g_Ykm1 double
        p_gamma_k double
        Xkp1f_est (:, 1, :) double
        Pkp1f (:, :, :) double
        Ykp1f_est (:, 1, :) double
        xkp1_est (:, 1) double
        Pkp1 double
        ykp1_est (:, 1) double
        xk_est (:, 1) double
        Pk double
        yk_est (:, 1) double
        i (1, 1) {mustBeInteger, mustBeNonnegative}
        type (1, 1) string
    end
    methods
        function obj = AbstractMKFObserverGPB(A,B,C,Ts,P0,Q,R,T,n_filt,label,x0, ...
                p_seq_g_Yk_init)

            % System dimensions
            [n, nu, ny] = check_dimensions(A{1}, B{1}, C{1});

            % Number of switching systems
            nj = numel(A);

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

            % Sequence index starts at 0 but is 1 thereafter.
            % TODO: Do we need this for comparibility?
            %obj.i = int16(0);

            % Fusion horizon length
            obj.f = 1;

            % Prior assumptions at initialization
            obj.p_seq_g_Yk_init = p_seq_g_Yk_init;

            % Check transition probability matrix
            % TODO: Is this necessary? If all possible hypotheses 
            % are not modelled (as is the case with some observers)
            % then perhaps T should not be whole?
            % assert(all(abs(sum(obj.T, 2) - 1) < 1e-15), "ValueError: T")
            % TODO: Could calculate beta parameter here, i.e. total 
            % probability measured?

            % Check all other system matrix dimensions have same 
            % input/output dimensions and number of states.
            for j = 2:nj
                [n_j, nu_j, ny_j] = check_dimensions(A{j}, B{j}, C{j});
                assert(isequal([n_j, nu_j, ny_j], [n, nu, ny]), ...
                    "ValueError: size of A, B, C, and D")
            end

            % Store useful variables
            obj.n = n;
            obj.nu = nu;
            obj.ny = ny;
            obj.nj = nj;
            obj.n_filt = n_filt;
            obj.type = "MKF_GPB";

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
            % At initialization time k = 0, xkp1_est, ykp1_est,
            % and Pkp1 represent prior estimates i.e. x_est(k|k-1),
            % y_est(k|k-1), and P(k|k-1).
            obj.xkp1_est = obj.x0;
            obj.ykp1_est = obj.C{1} * obj.xkp1_est;
            obj.Pkp1 = obj.P0;

            % At initialization time k = 0, x_est(k|k), y_est(k|k)
            % and P(k|k) have not been computed.
            obj.xk_est = nan(obj.n, 1);
            obj.yk_est = nan(obj.ny, 1);
            obj.Pk = nan(obj.n);

            % Initialise arrays to store predictions of all filters 
            obj.Xkp1f_est = repmat(obj.xkp1_est, 1, 1, obj.n_filt);
            obj.Pkp1f = repmat(obj.Pkp1, 1, 1, obj.n_filt);
            obj.Ykp1f_est = repmat(obj.ykp1_est, 1, 1, obj.n_filt);

        end
    end
end
