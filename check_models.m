function [nj, n, nu, ny, Ts] = check_models(models)
% [nj, n, nu, ny, Ts] = check_models(models)
% Checks all system models have same input/output 
% dimensions and number of states and returns the
% dimensions.
%
% Arguments:
%   models : (1, nj) cell array of structs
%       Each struct contains the parameters of a linear
%       dynamical system model. These include: A, B, 
%       and C for the system matrices, Q and R for the
%       state error covariance and output measurement 
%       noise covariance, and Ts for the sample period.
%

    % Number of models
    nj = numel(models);

    % Check and get dimensions of first model
    if isprop(models{1}, "D")
        [n, nu, ny] = check_dimensions(models{1}.A, models{1}.B, ...
                models{1}.C, models{1}.D);
    else
        [n, nu, ny] = check_dimensions(models{1}.A, models{1}.B, ...
                models{1}.C);
    end

    % Sample period
    Ts = models{1}.Ts;

    % Check other models have same dimensions.
    for j = 2:nj
        if isprop(models{j}, "D")
            [n_j, nu_j, ny_j] = check_dimensions(models{j}.A, ...
                models{j}.B, models{j}.C, models{j}.D);
        else
            [n_j, nu_j, ny_j] = check_dimensions(models{j}.A, ...
                models{j}.B, models{j}.C);
        end
        assert(isequal([n_j, nu_j, ny_j], [n, nu, ny]), ...
            "ValueError: size of A, B, and C")
        assert(models{j}.Ts == Ts, ...
            "ValueError: Ts")
    end

end