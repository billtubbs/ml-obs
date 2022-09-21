% Multi-model Kalman Filter class definition
%
% Class for simulating a multi-model Kalman filter for state
% estimation of a Markov jump linear system. This version
% implements a 'detection interval' used by certain sub-optimal
% algorithms to reduce the number of filters required. See
% MKFObserverSF.m for example.
%
% obs = MKFObserverDI(models,P0,seq,T,d,label,x0, ...
%     p_seq_g_Yk_init,reset)
%
% Arguments:
%	A, B, C : cell arrays containing discrete-time system
%       matrices for each switching system modelled.
%   Ts : sample period.
%   P0 : Initial covariance matrix of the state estimates
%       (same for each filter).
%   Q : Cell array of process noise covariance matrices for
%       each switching system.
%   R : Cell array of output measurement noise covariance
%       matrices for each switching system.
%   seq : Model indicator sequences for each filter (in rows).
%   T : Transition probabity matrix of the Markov switching
%       process.
%   d : Detection interval length in number of sample periods.
%   label : string name.
%   x0 : Initial state estimates (optional, default zeros).
%   r0 : (optional, default zeros)
%       Initial prior model indicator value at time k-1 
%       (zero-based, i.e. 0 is for first model).
%
% References:
%  -  Robertson, D. G., & Lee, J. H. (1998). A method for the
%     estimation of infrequent abrupt changes in nonlinear 
%     systems. Automatica, 34(2), 261-270.
%     https://doi.org/10.1016/S0005-1098(97)00192-1%
%

% TODO: This should inherit from MKFObserver to avoid duplication.

classdef MKFObserverDI < MKFObserverS
    properties (SetAccess = immutable)
        d (1, 1) double {mustBeInteger, mustBeNonnegative}
        f (1, 1) double {mustBeInteger, mustBeNonnegative}
        n_di (1, 1) double 
    end
    properties
        i2 int16
        i2_next int16
        idx_branch cell
        idx_modes cell
        idx_merge cell
    end
    methods
        function obj = MKFObserverDI(models,P0,seq,T,d,label,x0,r0, ...
                p_seq_g_Yk_init,reset)
            arguments
                models (1, :) cell
                P0 double
                seq (:, 1) cell
                T double
                d double
                label (1, 1) string = ""
                x0 = []
                r0 (:, 1) int16 {mustBeGreaterThan(r0, 0)} = 1
                p_seq_g_Yk_init = []
                reset logical = true
            end

            % Create super-class observer instance
            obj = obj@MKFObserverS(models,P0,seq,T,label, ...
                x0,r0,p_seq_g_Yk_init,false)

            % Number of detection intervals in horizon
            obj.n_di = size(cell2mat(obj.seq), 2);

            % Fusion horizon length in number of sample times
            obj.f = obj.n_di * d;

            % Generate the branch, mode, and merge index vectors 
            % - these govern the branching, mode transitions,
            % and merging of estimates at the end of each 
            % detection interval.
            [obj.idx_branch, obj.idx_modes, obj.idx_merge] = ...
                seq_fusion_indices(cell2mat(seq), obj.nj);

            % Add other useful variables
            obj.d = d;
            obj.type = "MKF_DI";

            if reset
                % Initialize variables
                obj.reset()
            end

        end
        function reset(obj)
        % obj.reset()
        % Initialize all variables to initial values specified
        % when observer object was created.
        %

            % Call reset method of super class object
            reset@MKFObserverS(obj)

            % Sequence index and counter for prob. updates
            % obj.i is the sequence index (1 <= i(1) <= obj.f)
            % obj.i2 is the counter for prob. updates (1 <= i(2) <= obj.d)
            obj.i2 = int16(0);
            obj.i2_next = int16(1);

        end
        function update(obj, yk, uk)
        % obj.update(yk, uk)
        % updates the estimates of the multi-model Kalman filter
        % and calculates the predictions of the states and output
        % at the next sample time.
        %
        % Arguments:
        %   obs : struct containing the multi-model Kalman filter
        %       variables (see function mkf_filter).
        %   yk : vector (ny, 1) of system output measurements
        %       at the current sample time.
        %   uk : vector (nu, 1) of system inputs at the current 
        %       sample time.
        %

            % Increment sequence index and update counter
            % obj.i is the sequence index (1 <= i <= obj.f)
            % obj.i2 is the counter for prob. updates (1 <= i2 <= obj.d)
            % Whenever obj.i2 exceeds obj.d (the spacing parameter), it is
            % reset to 1, and the sequence index obj.i is incremented.
            obj.i = obj.i_next;
            obj.i2 = obj.i2_next;
            obj.i_next = mod(obj.i - 1 + ...
                idivide(obj.i2, obj.d), obj.n_di) + 1;
            obj.i2_next = mod(obj.i2, obj.d) + 1;

            % Set model indicator values r(k) to the current
            % values from the filter's sequence and keep a copy 
            % of the previous values
            rkm1 = obj.rk;
            obj.rk = cell2mat(cellfun(@(x) x(:, obj.i), obj.seq, ...
                'UniformOutput', false));

            % TODO: Above can be simplified and do not need to do
            % Bayesian update every increment (but this should work)
            % Call update method of super class object
            update@MKFObserver(obj, yk, uk, obj.rk);

            % TODO: Do we need to implement branching and merging 
            % (should this observer inherit from MKFObserverSF instead?
%             % If at start of a detection interval, carry out
%             % filter branching and merging.
%             if obj.i2 == 1
%                 merge_idx = obj.merge_idx_seq(:, obj.i);
%                 counts = accumarray(merge_idx, 1);
%                 % Find hypothesis scheduled to be merged
%                 idxs_to_merge = find(counts > 1);
%                 for j = 1:numel(idxs_to_merge)
%                     %disp("Merging")
%                 end
%             end

        end
    end
end
