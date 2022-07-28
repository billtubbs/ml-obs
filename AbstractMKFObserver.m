% Multi-model observer class definition
%
% Abstract base class for the generalised pseudo-Bayes multi-
% model Kalman filters for state estimation of Markov 
% jump linear systems. The defining characteristic of
% these methods is that the number of independent
% Kalman filters matches the number of system models.
% I.e. obj.n_filt = obj.nj
%
% obs = AbstractMKFObserver(models,Ts,n_filt,label,x0,y0)
%
% Arguments:
%   models : cell array of structs defining a set of 
%       discrete-time linear system models. Each struct
%       must have the system matrices A, B, C, and the
%       covariance matrices Q and R as fields.
%   Ts : Sample period.
%   P0 : Initial covariance of the state estimation errors
%       (same for each filter).
%   n_filt : Number of independent filters to model.
%   label : string name.
%   x0 : Initial state estimates (optional, default zeros).
%   y0 : Initial output estimates (optional, default zeros).
%

classdef (Abstract) AbstractMKFObserver < matlab.mixin.Copyable
    properties (SetAccess = immutable)
        Ts (1, 1) double {mustBeNonnegative}
        n (1, 1) double {mustBeInteger, mustBeNonnegative}
        nu (1, 1) double {mustBeInteger, mustBeNonnegative}
        ny (1, 1) double {mustBeInteger, mustBeNonnegative}
        nj (1, 1) double {mustBeInteger, mustBeNonnegative}
        n_filt (1, 1) double {mustBeInteger, mustBeNonnegative}
    end
    properties
        models (1, :) cell
        label (1, 1) string
        x0 (:, 1) double
        P0 double
        y0 (:, 1) double
        xk_est (:, 1) double
        Pk double
        yk_est (:, 1) double
        type (1, 1) string
    end
    methods
        function obj = AbstractMKFObserver(models,Ts,n_filt,label,x0,y0)

            % Number of modes (switching systems)
            nj = numel(models);

            % Number of system inputs and outputs
            [n, nu, ny] = check_dimensions(models{1}.A, models{1}.B, ...
                models{1}.C);

            % Check all system models have same dimensions.
            for j = 2:nj
                if isfield(models{j}, "D")
                    % Check no direct transmission (D = 0)
                    assert(all(models{j}.D == 0, [1 2]), ...
                        "ValueError: D matrix")
                end
                [n2, nu2, ny2] = check_dimensions(models{j}.A, ...
                    models{j}.B, models{j}.C);
                assert(isequal([n2, nu2, ny2], [n, nu, ny]), ...
                    "ValueError: size of A, B, C matrices")
            end

            % Save properties
            obj.models = models;
            obj.Ts = Ts;
            obj.n_filt = n_filt;
            obj.label = label;
            obj.x0 = x0;
            obj.y0 = y0;
            obj.n = n;
            obj.nu = nu;
            obj.ny = ny;
            obj.nj = nj;
            obj.type = "MKF_GPB";

            % Initialize variables
            obj.reset();

        end
        function reset(obj)
        % obj.reset()
        % Initialize all variables to initial values specified
        % when observer object was created.
        %

            % At initialization time k = 0, x_est(k|k), y_est(k|k)
            % and P(k|k) are not defined.
            obj.xk_est = nan(obj.n, 1);
            obj.yk_est = nan(obj.ny, 1);
            obj.Pk = nan(obj.n);

        end
    end
end
