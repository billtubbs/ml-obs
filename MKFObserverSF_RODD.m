% Multi-model Kalman Filter class definition
%
% Object class for simulating a multi-model observer for 
% state estimation in the presence of randomly-occurring 
% deterministic disturbances (RODDs) as described in 
% Robertson et al. (1995, 1998).
%
% obs = MKFObserverSF_RODDMKFObserverSF_RODD(model,u_meas,P0,epsilon, ...
%     sigma_wp,Q0,R,f,m,d,label,x0)
%
% Arguments:
%   A, B, C : Matrices of the discrete time state-space
%       system representing the augmented system (including
%       disturbances and unmeasured inputs).
%   Ts : Sample period.
%   u_meas : binary vector indicating measured inputs.
%   P0 : Initial value of covariance matrix of the state
%       estimates.
%   epsilon : probability of a shock disturbance.
%   sigma_wp : standard deviation of shock disturbances.
%   Q0 : Process noise covariance matrix (n, n) with 
%        variances for each state on the diagonal. The  
%        values for states impacted by the unmeasured
%        input disturbances should be set to zero as the
%        appropriate variances will be added by the
%        algorithm during observer updates.
%   R : Output measurement noise covariance matrix (ny, ny).
%   f : Fusion horizon (length of disturbance sequences).
%   m : Maximum number of disturbances over fusion horizon.
%   d : Detection interval length in number of sample periods.
%   label : String name.
%   x0 : Initial state estimates (optional, default zeros)
%
% References:
%  -  Robertson, D. G., Kesavan, P., & Lee, J. H. (1995). 
%     Detection and estimation of randomly occurring 
%     deterministic disturbances. Proceedings of 1995 American
%     Control Conference - ACC'95, 6, 4453-4457. 
%     https://doi.org/10.1109/ACC.1995.532779
%  -  Robertson, D. G., & Lee, J. H. (1998). A method for the
%     estimation of infrequent abrupt changes in nonlinear 
%     systems. Automatica, 34(2), 261-270.
%     https://doi.org/10.1016/S0005-1098(97)00192-1%
%

classdef MKFObserverSF_RODD < MKFObserverDI
    properties (SetAccess = immutable)
        u_meas {mustBeNumericOrLogical}
        m double {mustBeInteger, mustBeNonnegative}
    end
    properties
        sys_model struct
        alpha double
        beta (1, 1) double
        p_seq double
        p_gamma double
        Q0 double
        R double
        epsilon double
        sigma_wp double
    end
    methods
        function obj = MKFObserverSF_RODD(model,u_meas,P0,epsilon, ...
                sigma_wp,Q0,R,f,m,d,label,x0)
            arguments
                model struct
                u_meas (:, 1) logical
                P0 double
                epsilon double
                sigma_wp double
                Q0 (:, :) double
                R (:, :) double
                f (1, 1) double {mustBeInteger}
                m (1, 1) double {mustBeInteger}
                d (1, 1) double {mustBeInteger}
                label (1, 1) string = ""
                x0 = []
            end

            % Number of states
            [n, nu, ny] = check_dimensions(model.A, model.B, model.C);
            if isprop(model, "D")
                assert(isequal(model.D, zeros(ny, nu)), ...
                    "ValueError: direct transmission (model.D)")
            end

            % Determine initial state values
            if isempty(x0)
                x0 = zeros(n, 1);  % default initial states
            else
                assert(isequal(size(x0), [n 1]), "ValueError: x0")
            end

            % Check size of initial process covariance matrix
            assert(isequal(size(Q0), [n n]), "ValueError: size(Q0)")

            % Observer model without disturbance noise input
            nw = sum(~u_meas);  % Number of input disturbances
            assert(nw > 0, "ValueError: u_meas");
            Bu = model.B(:, u_meas);
            Bw = model.B(:, ~u_meas);

            % Probability of at least one shock in a detection interval
            % (Detection interval is d sample periods in length).
            if d == 1
                alpha = epsilon;
            else
                alpha = (1 - (1 - epsilon).^d);
            end

            % Modified variances of shock signal over detection
            % interval (see (16) on p.264 of Robertson et al. 1998)
            var_wp = sigma_wp.^2;
            var_wp(~u_meas) = var_wp(~u_meas) ./ d;

            % Convert fusion horizon to number of detection intervals
            assert(rem(f, d) == 0, "ValueError: Fusion horizon and " ...
                + "detection interval not compatible")
            n_di = f / d;

            % Construct process noise covariance matrices and switching
            % sequences over the fusion horizon, and the prior 
            % probabilities of each sequence.
            [Q, p_rk, seq] = construct_Q_model_SF(Q0, Bw, alpha, ...
                var_wp, n_di, m, nw);

            % Number of models (each with a different hypothesis sequence)
            nj = numel(Q);

            % Transition probability matrix
            % Note that for RODD disturbances Pr(gamma(k)) is
            % assumed to be an independent random variable.
            T = repmat(p_rk', nj, 1);

            % Sequence probabilities Pr(Gamma(k))
            p_seq = prob_seq(seq, p_rk);

            % Tolerance parameter (total probability of defined sequences)
            beta = sum(p_seq);

            % Initial covariance matrix is the same for all filters
            %P0_init = repmat({P0}, 1, nh);

            % System models are all the same - only Q switches
            model_obs = model;
            model_obs.B = Bu;
            models = repmat({model_obs}, 1, nj);
            for i = 1:nj
                models{i}.Q = Q{i};
                models{i}.R = R;
            end

            % Number of hypotheses to be modelled
            nh = size(seq, 1);

            % Assume prior probabilities are equal (uniform)
            p_seq_g_Yk_init = ones(nh, 1) / nh;

            % Create MKF super-class observer instance
            obj = obj@MKFObserverDI(models,P0,seq,T,d,label,x0, ...
                p_seq_g_Yk_init,false)

            % Add additional variables used by RODD observer
            obj.sys_model = model;
            obj.u_meas = u_meas;
            obj.P0 = P0;
            obj.Q0 = Q0;
            obj.R = R;
            obj.epsilon = epsilon;
            obj.sigma_wp = sigma_wp;
            obj.m = m;
            obj.alpha = alpha;
            obj.beta = beta;
            obj.p_seq = p_seq;
            obj.p_gamma = p_rk;
            obj.type = "MKF_SF_RODD";

            % Initialize all variables
            obj.reset()

        end
    end
end
