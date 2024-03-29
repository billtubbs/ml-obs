function [Xk_m, Pk_m, Yk_m, p_seq_g_Yk_m] = merge_estimates(Xk, Pk, Yk, ...
    p_seq_g_Yk, rk, nj)
% [Xk_m, Pk_m, p_seq_g_Yk_m] = merge_estimates(Xk, Pk, ...
%     p_seq_g_Yk, rk, nj)
% Merge (i.e. fuse) a set of multiple-model Kalman filter
% estimates according to the modes specified in vector rk.
% Note that the sets of states, Xk, and covariance matrices,
% Pk, are both 3-dimensional arrays with the 3rd dimension
% represeting the system modes over which to merge.
%
% Example:
% >> Xk = cat(3, 1, 2, 2, 3);
% >> Pk = cat(3, 10, 11, 12, 13);
% >> Yk = cat(3, 2, 4, 4, 6);
% >> p_seq_g_Yk = [0.8 0.2 0.2 0.8]';
% >> rk = [1 1 2 2]'; nj = 2;
% >> [Xk_m,Pk_m,Yk_m,p_seq_g_Yk_m] = ...
%     merge_estimates(Xk,Pk,Yk,p_seq_g_Yk,rk, nj)
% 
% Xk_m(:,:,1) =
% 
%     1.2000
% 
% 
% Xk_m(:,:,2) =
% 
%     2.8000
% 
% 
% Pk_m(:,:,1) =
% 
%    10.3600
% 
% 
% Pk_m(:,:,2) =
% 
%    12.9600
% 
% 
% Yk_m(:,:,1) =
% 
%     2.4000
% 
% 
% Yk_m(:,:,2) =
% 
%     5.6000
% 
% 
% p_seq_g_Yk_m =
% 
%      1
%      1
%

    % Number of states
    n = size(Xk, 1);

    % Number of Outputs
    ny = size(Yk, 1);

    % Number of hypothesis
    nh = size(p_seq_g_Yk, 1);

    % Weighting to apply to hypotheses
    weights = reshape(p_seq_g_Yk, 1, 1, nh);

    if nargin == 4
        % If p_seq_g_Yk, rk, and nj are not provided, do
        % a complete merge to one set of estimates

        % Merge state and output estimates
        Xk_m = sum(weights .* Xk, 3);
        Yk_m = sum(weights .* Yk, 3);

        % Merge state estimation error covariance
        X_diffs = Xk_m - Xk;
        Pk_m = sum(weights .* (Pk + ...
            pagemtimes(X_diffs, pagetranspose(X_diffs))), 3);

        % Merged probabilities (sum to 1)
        p_seq_g_Yk_m = 1;

    else

        % Partial merge based on index vector rk

        % Create mask to allow vectorized merging calculation
        mask = false(1, 1, nh, nj);
        mask(sub2ind(size(mask), ones(nh, 1), ones(nh, 1), (1:nh)', rk)) = true;

        % Normalisation constants
        C = sum(weights .* mask, 3);

        % Merge state and output estimates
        Xk_m = reshape(sum(weights .* Xk .* mask, 3) ./ C, n, 1, []);
        Yk_m = reshape(sum(weights .* Yk .* mask, 3) ./ C, ny, 1, []);

        % Merge state estimation error covariance
        X_diffs = Xk_m(:, :, rk) - Xk;
        to_add = pagemtimes(X_diffs, pagetranspose(X_diffs));
        Pk_m = reshape(sum(weights .* (Pk + to_add) .* mask, 3) ./ C, n, n, []);

        % Merged probabilities
        p_seq_g_Yk_m = reshape(C, [nj 1]);

    end

end