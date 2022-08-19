function [xk_est,yk_est,Pk,p_seq_g_Yk] = merge_estimates(xk_est, Pk, yk_est, ...
    p_seq_g_Yk, gamma_k, nj)
% [xk_est,yk_est,Pk,p_seq_g_Yk] = merge_estimates(xk_est, Pk, yk_est, p_seq_g_Yk, ...
%     gamma_k, nj)
% Merge (i.e. fuse) a set of multiple-model Kalman filter
% estimates according to the modes specified in gamma_k.
%
% Example:
% >> xk_est = [1 2 2 3]';
% >> p_seq_g_Yk = [0.8 0.2 0.2 0.8]';
% >> gamma_k = [0 0 1 1]';
% >> merge_estimates(xk_est, p_seq_g_Yk, gamma_k, 2)
% 
% ans =
% 
%     1.2000
%     2.8000
%
    n_filt = size(gamma_k, 1);
    n = size(xk_est, 1);
    ny = size(yk_est, 1);
    mask = false(1, n_filt, nj);
    mask(sub2ind(size(mask), ones(n_filt, 1), (1:n_filt)', gamma_k + 1)) = true;
    xk_est_m = reshape(sum((xk_est .* p_seq_g_Yk') .* mask, 2), [n nj]);
    yk_est_m = reshape(sum((yk_est .* p_seq_g_Yk') .* mask, 2), [ny nj]);
    Pk_m = zeros(n, n, nj);
    for f = 1:n_filt
        m_ind = gamma_k(f) + 1;
        %xk_est(:, m_ind) = xk_est(:, m_ind) + updx(:, f) .* p_seq_g_Yk(f);
        %yk_est(:, m_ind) = yk_est(:, m_ind) + updy(:, f) .* p_seq_g_Yk(f);
        x_diff = xk_est_m(:, m_ind) - xk_est(:, f);
        Pk_m(:,:,m_ind) = Pk_m(:,:,m_ind) ...
            + p_seq_g_Yk(f) * (Pk(:, :, f) + x_diff * x_diff');
    end
    xk_est = xk_est_m;
    yk_est = yk_est_m;
    Pk = Pk_m;

end