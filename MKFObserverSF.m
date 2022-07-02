% Multi-model Kalman Filter class definition
%
% Object class for simulating a multi-model observer for 
% state estimation in the presence of randomly-occurring 
% deterministic disturbances (RODDs) as described in 
% Robertson et al. (1995, 1998).
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

classdef MKFObserverSF < MKFObserverDI
    properties (SetAccess = immutable)
        u_meas {mustBeNumericOrLogical}
        m double {mustBeInteger, mustBeNonnegative}
    end
    properties
        alpha double
        beta (1, 1) double
        p_seq double
        p_gamma double
        Q0 {mustBeNumeric}
        epsilon double
        sigma_wp double
    end
    methods
        function obj = MKFObserverSF(A,B,C,D,Ts,u_meas,P0,epsilon, ...
                sigma_wp,Q0,R,f,m,d,label,x0)
        % obs = MKFObserverSF(A,B,C,D,Ts,u_meas,P0,epsilon, ...
        %     sigma_wp,Q0,R,f,m,d,label,x0)
        %
        % Arguments:
        %   A, B, C, D : matrices of the discrete time state-space
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
        %   x0 : intial state estimates (optional).
        %

            % Number of states
            n = check_dimensions(A, B, C, D);

            % Check size of initial process covariance matrix
            assert(isequal(size(Q0), [n n]), "ValueError: size(Q0)")

            % Initial state estimates
            if nargin < 16
                x0 = zeros(n,1);
            end

            % Check size of process covariance default matrix
            assert(isequal(size(Q0), [n n]), "ValueError: size(Q0)")

            % Observer model without disturbance noise input
            nw = sum(~u_meas);  % Number of input disturbances
            assert(nw > 0, "ValueError: u_meas");
            Bu = B(:, u_meas);
            Bw = B(:, ~u_meas);
            Du = D(:, u_meas);

            % Probability of at least one shock in a detection interval
            % (Detection interval is d sample periods in length).
            if d == 1
                alpha = epsilon;
            else
                alpha = (1 - (1 - epsilon).^d);
            end

            % Modified variances of shock signal over detection
            % interval (see (16) on p.264 of Robertson et al. 1998)
            var_wp = sigma_wp.^2 ./ d;

            % Convert fusion horizon to number of detection intervals
            assert(rem(f, d) == 0, "ValueError: Fusion horizon and " ...
                + "detection interval not compatible")
            n_di = f / d;

            % Construct process noise covariance matrices and switching
            % sequences over the fusion horizon, and the prior 
            % probabilities of each sequence.
            [Q, p_gamma, seq] = construct_Q_model_SF(Q0, Bw, alpha, ...
                var_wp, n_di, m, nw);

            % Number of models (each with a different hypothesis sequence)
            nj = numel(Q);

            % Transition probability matrix
            % Note that for RODD disturbances Pr(gamma(k)) is
            % assumed to be an independent random variable.
            T = repmat(p_gamma', nj, 1);

            % Sequence probabilities Pr(Gamma(k))
            p_seq = prob_Gammas(seq, p_gamma);

            % Tolerance parameter (total probability of defined sequences)
            beta = sum(p_seq);

            % System model doesn't change
            A = repmat({A}, 1, nj);
            Bu = repmat({Bu}, 1, nj);
            C = repmat({C}, 1, nj);
            Du = repmat({Du}, 1, nj);
            R = repmat({R}, 1, nj);

            % Initial covariance matrix is the same for all filters
            %P0_init = repmat({P0}, 1, n_filt);

            % Create MKF super-class observer instance
            obj = obj@MKFObserverDI(A,Bu,C,Du,Ts,P0,Q,R,seq,T,d,label,x0);

            % Add additional variables used by RODD observer
            obj.u_meas = u_meas;
            obj.P0 = P0;
            obj.Q0 = Q0;
            obj.epsilon = epsilon;
            obj.sigma_wp = sigma_wp;
            obj.m = int16(m);
            obj.alpha = alpha;
            obj.beta = beta;
            obj.p_seq = p_seq;
            obj.p_gamma = p_gamma;
            obj.type = "MKF_SF";

        end
    end
end
