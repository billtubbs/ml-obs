% Kalman Filter class definition
%
% obs = SKFObserverS(models,P0,seq,label,x0,r0)
% Class for simulating a dynamic Kalman filter
% (i.e. with time-varying gain and estimation error
% covariance) for state estimation of a jump linear
% system where the system model switches between
% a finite set of models each time step.
%
% This version differs from SKFObserver in that it
% has a pre-determined sequence of system modes
% (i.e. transitions) that it models.
%
% This is the filtering form of the observer which 
% produces posterior estimates of the states and
% outputs at the current time instant given the data
% at the current time:
%
%   x_hat(k|k) : estimate of states at time k
%   y_hat(k|k) : estimate of outputs at time k
%
% Prior estimates of the states at the next time
% instant given the data at the current time are also
% calculated:
%
%   x_hat(k+1|k) : estimate of states at time k + 1
%
% The system model is defined as:
%
%   x(k+1) = A(k) x(k) + B(k) u(k) + w(k)
%     y(k) = C(k) x(k) + v(k)
%
% Note: there is no direct transmission (D = 0).
%
% Arguments:
%   models : (1, nj) cell array of structs
%       Each struct contains the parameters of a linear
%       model of the system dynamics. These include: A, B, 
%       and C for the system matrices, Q and R for the
%       state error covariance and output measurement 
%       noise covariance, and Ts for the sample period.
%   P0 : (n, n) double
%       Initial covariance matrix of the state estimates
%       (same for each filter).
%   seq : (1, nf) double
%       Switching sequence of system (i.e. mode transitions)
%       to be modelled. These should be integers in the 
%       range {1, ..., nj} which indicate the system mode
%       at time k = 0, 1, ..., nf-1.
%   label : string
%       Arbitrary name to identify the observer.
%   x0 : (n, 1) double (optional, default zeros)
%       Initial state estimates.
%   r0 : (nh, 1) integer (optional, default ones)
%       Prior system mode at time k = -1.
%

classdef SKFObserverS < SKFObserver
    properties
        seq (1, :) double
        nf (1, 1) double {mustBeInteger}
        i (1, 1) {mustBeInteger, mustBeNonnegative}
        i_next (1, 1) {mustBeInteger, mustBeNonnegative}
    end
    methods
        function obj = SKFObserverS(models,P0,seq,label,x0,r0)
            arguments
                models (1, :) cell
                P0 double
                seq (1, :) double
                label (1, 1) string = ""
                x0 = []
                r0 = 1  % default system mode at k = 0
            end

            % Create super-class observer instance
            obj = obj@SKFObserver(models,P0,label,x0,r0);

            % Store parameters
            obj.seq = seq;
            obj.nf = size(seq, 2);
            obj.type = "SKF-S";

            % Initialize variables
            obj.reset()

        end
        function reset(obj)
        % obj.reset()
        % Initialize all variables to initial values specified
        % when observer object was created.
        %

            % Call reset method of super class object
            reset@SKFObserver(obj);

            % Reset sequence index
            obj.i = int16(0);
            obj.i_next = int16(1);

        end
        function update(obj, yk, uk)
        % obj.update(yk, uk, rk) updates the gain and covariance 
        % matrix of the Kalman filter and calculates the estimates 
        % of the states and output at the next sample time given
        % the current measurement, inputs, and system mode.
        %
        % Arguments:
        %   yk : vector, size (ny, 1)
        %       System output measurements at current time k.
        %   uk : vector, size (nu, 1)
        %       System inputs at current time k.
        %   rk : integer
        %       System mode at current time k.
        %

            % Increment sequence index (at end of sequence it 
            % wraps to beginnning)
            obj.i = obj.i_next;
            obj.i_next = mod(obj.i, obj.nf) + 1;

            % Get current system mode from sequence
            rk = obj.seq(:, obj.i);

            % Call reset method of super class object
            update@SKFObserver(obj, yk, uk, rk);

        end
    end
end