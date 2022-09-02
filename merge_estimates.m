function [xk_est_m,yk_est_m,Pk_m,p_seq_g_Yk_m] = merge_estimates( ...
    xk_est,Pk,yk_est,p_seq_g_Yk,gamma_k,nj)
% [xk_est_m,yk_est_m,Pk_m,p_seq_g_Yk_m] = merge_estimates(xk_est,Pk, ...
%     yk_est,p_seq_g_Yk,gamma_k,nj)
% Merge (i.e. fuse) a set of multiple-model Kalman filter
% estimates according to the modes specified in gamma_k.
%
% Example:
% >> xk_est = [1 2 2 3];
% >> Pk = cat(3, 10, 11, 12, 13);
% >> yk_est = [2 4 4 6];
% >> p_seq_g_Yk = [0.8 0.2 0.2 0.8]';
% >> gamma_k = [0 0 1 1]';
% >> [xk_est_m,yk_est_m,Pk_m,p_seq_g_Yk_m] = ...
%     merge_estimates(xk_est,Pk,yk_est,p_seq_g_Yk,gamma_k, nj)
% 
% xk_est_m =
% 
%     1.2000    2.8000
% 
% 
% yk_est_m =
% 
%     2.4000    5.6000
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
% p_seq_g_Yk_m =
% 
%      1
%      1
%

    n_filt = size(gamma_k, 1);
    n = size(xk_est, 1);
    ny = size(yk_est, 1);
    mask = false(1, n_filt, nj);
    mask(sub2ind(size(mask), ones(n_filt, 1), (1:n_filt)', gamma_k + 1)) = true;

    % Normalisation constants
    C = sum(p_seq_g_Yk' .* mask, 2);

    % Merge state and output estimates
    xk_est_m = reshape(sum((xk_est .* p_seq_g_Yk') .* mask, 2) ./ C, [n nj]);
    yk_est_m = reshape(sum((yk_est .* p_seq_g_Yk') .* mask, 2) ./ C, [ny nj]);

    % Merge state error covariance
    Pk_m = zeros(n, n, nj);
    for f = 1:n_filt
        m_ind = gamma_k(f) + 1;
        x_diff = xk_est_m(:,m_ind) - xk_est(:,f);
        Pk_m(:,:,m_ind) = Pk_m(:,:,m_ind) ...
            + p_seq_g_Yk(f) * (Pk(:, :, f) + x_diff * x_diff');
    end
    Pk_m = Pk_m ./ C;

    % Merged probabilities
    p_seq_g_Yk_m = reshape(C, [nj 1]);

end