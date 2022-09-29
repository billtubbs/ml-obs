% Multi-model Kalman Filter class definition
%
% Object class for simulating a multi-model observer for 
% state estimation in the presence of randomly-occurring 
% deterministic disturbances (RODDs) as described in 
% Robertson et al. (1995). This version is slightly 
% different to that described in Robertson et al. (1998).
%
% obs = MKFObserverSF95(A,B,C,Ts,u_meas,P0,epsilon, ...
%               sigma_wp,Q0,R,f,m,d,label,x0)
%
% Arguments:
%   A, B, C : matrices of the discrete time state-space
%       system representing the augmented system (including
%       disturbances and unmeasured inputs).
%   Ts : sample period.
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
%   R : output measurement noise covariance matrix (ny, ny).
%   f : fusion horizon (length of disturbance sequences).
%   m : maximum number of disturbances over fusion horizon.
%   d : detection interval length in number of sample periods.
%   label : string name.
%   x0 : Initial state estimates (optional, default zeros).
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

classdef MKFObserverSF_RODD95 < MKFObserverS
    properties (SetAccess = immutable)
        u_meas {mustBeNumericOrLogical}
        m double {mustBeInteger, mustBeNonnegative}
        d (1, 1) double {mustBeInteger, mustBeNonnegative}
    end
    properties
        sys_model struct
        beta (1, 1) double
        p_seq double
        p_rk double
        Q0 {mustBeNumeric}
        R double
        epsilon double
        sigma_wp double
    end
    methods
        function obj = MKFObserverSF_RODD95(model,u_meas,P0,epsilon, ...
                sigma_wp,Q0,R,f,m,d,label,x0,r0)
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
                r0 (:, 1) int16 {mustBeGreaterThan(r0, 0)} = 1
            end

            % Get model dimensions and sample period
            [n, nu, ny, Ts, direct] = check_model(model);
            if direct
                assert(isequal(model.D, zeros(ny, nu)), ...
                    "ValueError: direct transmission (model.D)")
            end

            % Determine initial state values
            if isempty(x0)
                x0 = zeros(n, 1);  % default initial states
            else
                assert(isequal(size(x0), [n 1]), "ValueError: x0")
            end

            % Check size of process covariance default matrix
            assert(isequal(size(Q0), [n n]), "ValueError: size(Q0)")

            % Observer model without disturbance noise input
            nw = sum(~u_meas);  % Number of input disturbances
            assert(nw > 0, "ValueError: u_meas");
            Bu = model.B(:, u_meas);
            Bw = model.B(:, ~u_meas);

            % Convert fusion horizon to number of detection intervals
            assert(rem(f, d) == 0, "ValueError: Fusion horizon and " ...
                + "detection interval not compatible")
            n_di = f / d;

            % Construct process noise covariance matrices and switching
            % sequences over the fusion horizon, and the prior 
            % probabilities of each sequence.
            [Q, p_rk, S] = construct_Q_model_SF95(Q0, Bw, epsilon, ...
                sigma_wp, n_di, m, nw);

            % Number of models (each with a different hypothesis sequence)
            nj = numel(Q);

            % Number of hypotheses (filters) to be modelled
            nh = size(S, 1);

            % Expand sequences by inserting zeros between times
            % when shocks occur.
            seq = cell(nh, 1);
            for i = 1:nh
                seq{i} = int16(ones(size(S{i}, 1), f));
                % Add shock indications at start of each detection
                % interval
                seq{i}(:, 1:d:f) = S{i};
                % Alternatively, at end of each detection interval
                %seq{i}(:, d:d:f) = S{i};
            end

            % Transition probability matrix
            % Note that for RODD disturbances Pr(gamma(k)) is
            % assumed to be an independent random variable.
            T = repmat(p_rk', nj, 1);

            % Sequence probabilities Pr(Gamma(k))
            p_seq = prob_seq(seq, p_rk);

            % Tolerance parameter (total probability of defined sequences)
            beta = sum(p_seq);

            % System models are all the same - only Q switches
            model_obs = model;
            model_obs.B = Bu;
            models = repmat({model_obs}, 1, nj);
            for i = 1:nj
                models{i}.Q = Q{i};
                models{i}.R = R;
            end

            % Only the Q parameter is different
            models{1}.Q = Q{1};
            models{2}.Q = Q{2};

            % Create MKF super-class observer instance
            obj = obj@MKFObserverS(models,P0,seq,T,label,x0);

            % Add additional variables used by RODD observer
            obj.sys_model = model;
            obj.u_meas = u_meas;
            obj.P0 = P0;
            obj.Q0 = Q0;
            obj.R = R;
            obj.epsilon = epsilon;
            obj.sigma_wp = sigma_wp;
            obj.m = m;
            obj.d = d;
            obj.beta = beta;
            obj.p_seq = p_seq;
            obj.p_rk = p_rk;
            obj.type = "MKF_SF_RODD95";

        end
    end
end
