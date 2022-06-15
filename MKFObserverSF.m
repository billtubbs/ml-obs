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

classdef MKFObserverSF < MKFObserver
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
        % obs = mkf_observer_RODD(A,B,C,D,Ts,u_meas,P0,epsilon, ...
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

            % Initial state estimates
            if nargin < 16
                x0 = zeros(n,1);
            end

            % Check size of process covariance default matrix
            assert(isequal(size(Q0), [n n]), "ValueError: size(Q0)")

            % Observer model without disturbance noise input
            Bu = B(:, u_meas);
            Bw = B(:, ~u_meas);
            Du = D(:, u_meas);
            nw = int16(sum(~u_meas));  % Number of input disturbances
            assert(nw > 0, "ValueError: u_meas");

            % Number of filters needed
            % TODO: arg f here is not actually the fusion horizon 
            % (which is f*d).  Should maybe use f/d her and assert no
            % remainder.
            n_filt = n_filters(m, f, nw);

            % Generate indicator sequences
            seq = combinations_lte(f*nw, m);

            % Probability of shock over a detection interval
            % (Detection interval is d sample periods in length).
            alpha = (1 - (1 - epsilon).^d);

            % Probabilities of no-shock / shock over detection interval
            % (this is named delta in Robertson et al. 1998)
            p_gamma = [1-alpha'; alpha'];

            if nw == 1

                % Number of models (with different Q matrices)
                nj = 2;

                % Generate required Q matrices.
                Q = cell(1, nj);
                for i = 1:nj
                    var_x = diag(Q0);
                    % Modified variances of shock signal over detection
                    % interval (see (16) on p.264 of Robertson et al. 1998)
                    var_x = var_x + Bw * sigma_wp(:, i).^2' ./ d;
                    Q{i} = diag(var_x);
                end

            elseif nw > 1

                % Note: In the case of more than one input disturbance
                % there may be multiple combinations of disturbances
                % occuring simultaneously. To simulate these, construct
                % a different Q matrix for each possible combination.

                % Reshape sequences into matrices with a row for each
                % input disturbance sequence
                seq = cellfun(@(x) reshape(x, nw, []), seq, ...
                    'UniformOutput', false);

                % Find all unique combinations of simultaneous shocks
                [Z,~,ic] = unique(cell2mat(seq')', 'sorted', 'rows');

                % Number of Q matrices needed
                nj = size(Z, 1);

                % Rearrange as one sequence for each filter and convert
                % back to cell array
                seq = reshape((ic - 1)', [], n_filt)'; 
                seq = mat2cell(seq, ones(n_filt, 1), f);

                % Generate required Q matrices
                Q = cell(1, nj);
                for i = 1:nj
                    ind = Z(i, :) + 1;
                    var_x = diag(Q0);
                    % Modified variances of shock signal over detection
                    % interval (see (16) on p.264 of Robertson et al. 1998)
                    idx = sub2ind(size(sigma_wp), 1:nw, ind);
                    var_x = var_x + Bw * sigma_wp(idx).^2' ./ d;
                    %var_x(~u_meas) = var_x(~u_meas) .* sigma_wp(idx).^2' ./ d;
                    Q{i} = diag(var_x);
                end

                % Modify indicator value probabilities for
                % combinations of shocks
                p_gamma = prod(prob_gamma(Z', p_gamma), 1)';

                % Normalize so that sum(Pr(gamma(k))) = 1
                % TODO: Is this the right thing to do for sub-optimal approach?
                p_gamma = p_gamma ./ sum(p_gamma);

            else

                error("Value error: no unmeasured inputs")

            end

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
            obj = obj@MKFObserver(A,Bu,C,Du,Ts,P0,Q,R,seq,T,d,label,x0);

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
        % TODO: Is MKFObserver reset method okay?
%         function reset(obj)
% 
%             disp("Not yet implemented")
%             reset@MKFObserver(obj)
% 
%         end

    end

end