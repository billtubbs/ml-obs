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


classdef MKFObserverGPB1 < AbstractMKFObserverGPB
    methods
        function obj = MKFObserverGPB1(A,B,C,Ts,P0,Q,R,T,label,x0, ...
                p_seq_g_Yk_init)

            % System dimensions
            n = check_dimensions(A{1}, B{1}, C{1});

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
            if nargin < 9
                label = "GPB1";
            end

            obj = obj@AbstractMKFObserverGPB(A,B,C,Ts,P0,Q,R,T,label,x0, ...
                p_seq_g_Yk_init);

            % Mode sequences (length 1 for GPB1)
            % TODO: Do we need this For compatibility with other 
            % MKF observers?
            %obj.seq = num2cell(obj.gamma_k);

            % Change type code
            obj.type = "MKF_GPB1";

            % Initialize all variables
            obj.reset()

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

            % Update state and output estimates based on current
            % measurement
            [obj.xk_est, obj.yk_est, obj.Pk, obj.p_seq_g_Yk] = ...
                GPB1_update(obj.A, obj.B, obj.C, obj.Q, obj.R, obj.T, ...
                    obj.Xkp1f_est, obj.Ykp1f_est, obj.Pkp1f, yk, ...
                    obj.p_seq_g_Yk);

            % Calculate predictions of each filter
            for j = 1:obj.n_filt
                [obj.Xkp1f_est(:,:,j), obj.Ykp1f_est(:,:,j), ...
                 obj.Pkp1f(:,:,j)] = kalman_predict_f(obj.A{j}, ...
                    obj.B{j}, obj.C{j}, obj.Q{j}, obj.xk_est, obj.Pk, uk);
            end

        end

    end

end
