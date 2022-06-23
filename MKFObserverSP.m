% Multi-model Kalman Filter class definition
%
% Object class for simulating a multi-model observer
% using the adaptive forgetting through multiple models 
% (AFMM) algorithm for state estimation in the presence of 
% infrequently-occurring deterministic disturbances, as 
% described in Eriksson and Isaksson (1996).
%
% Uses a sequence pruning method described in Eriksson and
% Isaksson (1996):
% - At the updates, let only the most probable sequence,
%   i.e. with the largest weight of all the sequences
%   split into 2 branches.
% - After the update is completed, remove the sequence
%   with the smallest weight.
% - Renormalize the remaining weights to have unit sum.
%
% Restriction to above rules:
% - Do not cut branches immediately after they are born.
%   Let there be a certain minimum life length for all 
%   branches.
%
% NOTE:
% - The adaptive forgetting component of the AFMM
%   (Andersson, 1985) is not yet implemented.
%
% References:
%  - Eriksson, P.-G., & Isaksson, A. J. (1996). Classification
%    of Infrequent Disturbances. IFAC Proceedings Volumes, 29(1), 
%     6614-6619. https://doi.org/10.1016/S1474-6670(17)58744-3
%  - Andersson, P. (1985). Adaptive forgetting in recursive
%    identification through multiple models?. International
%    Journal of Control, 42(5), 1175?1193. 
%    https://doi.org/10.1080/00207178508933420
%

classdef MKFObserverSP < MKFObserver
    properties (SetAccess = immutable)
        u_meas {mustBeNumericOrLogical}
        n_min double {mustBeInteger, mustBeNonnegative}
    end
    properties
        p_seq double
        p_gamma double
        Q0 {mustBeNumeric}
        epsilon double
        sigma_wp double
        n_main
        n_hold
        f_main
        f_hold
    end
    methods
        function obj = MKFObserverSP(A,B,C,D,Ts,u_meas,P0,epsilon, ...
            sigma_wp,Q0,R,n_filt,f,n_min,label,x0)
        % obs = MKFObserverSP(A,B,C,D,Ts,u_meas,P0,epsilon, ...
        %     sigma_wp,Q0,R,n_filt,f,n_min,label,x0)
        %
        % Creates a struct for simulating a multi-model observer
        % using the adaptive forgetting through multiple models 
        % (AFMM) method for state estimation in the presence of 
        % infrequently-occurring deterministic disturbances, as 
        % described in Eriksson and Isaksson (1996).
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
        %   n_filt : number of models (Kalman filters) to utilise.
        %   f : length of disturbance sequences to record.
        %   n_min : minimum life of cloned filters in number of
        %       sample periods.
        %   label : string name.
        %   x0 : intial state estimates (optional).
        %

            % TODO: Expose spacing parameter d as an argument
            % Number of states
            n = check_dimensions(A, B, C, D);

            % Check size of initial process covariance matrix
            assert(isequal(size(Q0), [n n]), "ValueError: size(Q0)")

            % Detection interval length in number of sample periods.
            d = 1;  % TODO: Make this a specified parameter.

            % Initial state estimates
            if nargin == 15
                x0 = zeros(n,1);
            end

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

            % Construct process noise covariance matrices for each 
            % possible input disturbance (returns a cell array)
            [Q, p_gamma] = construct_Q_model_SP(Q0, Bw, alpha, var_wp, nw);

            % Number of models (each with a different hypothesis sequence)
            nj = numel(Q);

            % Transition probability matrix
            % Note that for RODD disturbances Pr(gamma(k)) is
            % assumed to be an independent random variable.
            T = repmat(p_gamma', nj, 1);

            % Initialize indicator sequences
            seq = mat2cell(int16(zeros(n_filt, f)), int16(ones(1, n_filt)), f);

            % Define filter groups ('main', 'holding' and 'unused')
            n_min = int16(n_min);
            assert(n_min > 0)
            assert(n_filt > 0)
            n_hold = nw*n_min;
            n_main = n_filt - n_hold;
            
            % Check there are enough filters in total to accommodate
            % those in the holding group + at least one in main group
            assert(n_main >= nw, "ValueError: n_filt is too low.")

            % Filter indices
            f_main = int16(1:n_main);
            f_hold = int16(n_main+1:n_main+n_hold);

            % TODO: Remove this
            assert(isequal(sort(unique([f_main f_hold])), 1:(n_main+n_hold)))

            % System model doesn't change over time
            A = repmat({A}, 1, nj);
            Bu = repmat({Bu}, 1, nj);
            C = repmat({C}, 1, nj);
            Du = repmat({Du}, 1, nj);
            R = repmat({R}, 1, nj);

            % Initial covariance matrix is the same for all filters
            %P0_init = repmat({P0}, 1, n_filt);

            % Create MKF super-class observer instance
            obj = obj@MKFObserver(A,Bu,C,Du,Ts,P0,Q,R,seq,T,d,label,x0);

            % Sequence pruning algorithm initialization
            % Assign all probability to first filter
            obj.p_seq_g_Yk = [1; zeros(obj.n_filt-1, 1)];

            % Set estimate covariances to high values for the
            % rest of the filters
            for i = 2:obj.n_filt
                obj.filters{i}.P = 1e10*eye(obj.n);
            end

            % Add additional variables used by AFMM observer
            obj.u_meas = u_meas;
            obj.n_min = n_min;
            obj.n_main = n_main;
            obj.n_hold = n_hold;
            obj.f_main = f_main;
            obj.f_hold = f_hold;
            obj.P0 = P0;
            obj.Q0 = Q0;
            obj.epsilon = epsilon;
            obj.sigma_wp = sigma_wp;
            obj.p_gamma = p_gamma;
            obj.type = "MKF_SP";
        end
        function update(obj, yk, uk)
        % obs.update(yk, uk)
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

            % Implementation of filter pruning algorithm
            %
            % - As described in Eriksson and Isaksson (1996):
            %
            %   "Do not cut branches immediately after they are 
            %    born. Let there be a certain minimum life length 
            %    for all branches."
            %
            % To implement this rule, filters are organized into
            % two groups:
            % 
            %  1. Main group : obj.filters{f} for f in obj.f_hold
            %     Longer-surviving filters which can be killed off.
            %  2. Holding group : obj.filters{f} for f in obj.f_main
            %     Recently created filters which cannot be killed 
            %     off and are held for n_min time periods before
            %     being transferred to the main group.
            %
            % Note that for n_dist input disturbances, n_dist new
            % sequences are created each time step, so the
            % holding group needs to be n_dist*obj.n_min in size.

            % Index of current most likely sequence
            [~, f_max] = max(obj.p_seq_g_Yk);

            % Set next sequence value to 0 (no shock) for all
            % sequences
            for f = 1:obj.n_filt
                obj.seq{f}(:, obj.i_next(1)) = 0;
            end

            % Number of disturbances
            nw = size(obj.epsilon, 1);

            % Left-shift all filters in holding group. This causes
            % the last nw values to 'roll-over' to the left of f_hold.
            % e.g. circshift([1 2 3], 1) -> [3 1 2]
            obj.f_hold = circshift(obj.f_hold, -nw);

            % Filters to be moved out of holding group
            f_move = obj.f_hold(end-nw+1:end);

            % Rank sequences in main group according to 
            % conditional probabilities
            [~, i_rank] = sort(obj.p_seq_g_Yk(obj.f_main), 'descend');
            obj.f_main = obj.f_main(i_rank);

            % Select those with lowest probability for pruning
            f_to_prune = obj.f_main(end-nw+1:end);

            % Replace pruned sequences with those from holding
            % group
            obj.f_main(end-nw+1:end) = f_move;

            % Make clone(s) of most probable sequence and fitler
            % put new filter(s) at start of holding group
            obj.f_hold(end-nw+1:end) = f_to_prune;
            for i = 1:nw
                label = obj.filters{f_to_prune(i)}.label;
                obj.filters{f_to_prune(i)} = obj.filters{f_max}.copy();
                obj.filters{f_to_prune(i)}.label = label;  % keep label
                obj.p_seq_g_Yk(f_to_prune(i)) = obj.p_seq_g_Yk(f_max);
                obj.p_gamma_k(f_to_prune(i)) = obj.p_gamma_k(f_max);
                obj.p_seq_g_Ykm1(f_to_prune(i)) = obj.p_seq_g_Ykm1(f_max);
                obj.seq{f_to_prune(i)} = obj.seq{f_max};
                % Set next sequence value to index of the random input
                obj.seq{f_to_prune(i)}(:, obj.i_next(1)) = i;
            end

            % TODO: Add online noise variance estimation with
            % forgetting factor as described in Anderson (1995)
            % equation (4.3).

            % Run MKF super-class updates and probability calcs
            update@MKFObserver(obj, yk, uk);

        end
    end
end
