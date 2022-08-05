function [xk_est,yk_est,Pk,p_seq_g_Yk] = merge_estimates(xk_est, Pk, yk_est, p_seq_g_Yk, ...
    gamma_k, nj)
% [xk_est,yk_est,Pk,p_seq_g_Yk] = merge_estimates(xk_est, Pk, yk_est, p_seq_g_Yk, ...
%     gamma_k, nj)
% Merges (i.e. fuses) a set of multi-model estimates according
% to the modes specified in gamma_k.
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
    mask = false(n_filt, nj);
    mask(sub2ind(size(mask), (1:n_filt)', gamma_k + 1)) = true;
    xk_est = sum((xk_est .* p_seq_g_Yk) .* mask, 1)';
    yk_est = sum((yk_est .* p_seq_g_Yk) .* mask, 1)';
    Pk = zeros(n, n, nj);
    for f = 1:n_filt
        m_ind = gamma_k(f) + 1;
        xk_est(:, m_ind) = xk_est(:, m_ind) + updx(:, f) .* p_seq_g_Yk(f);
        yk_est(:, m_ind) = yk_est(:, m_ind) + updy(:, f) .* p_seq_g_Yk(f);
        Pk(:,:,m_ind) = Pk(:,:,m_ind) + p_seq_g_Yk(f) * (updP(:, :, f) + ...
            (updx(:, f) - xk_est) * (updx(:, f) - xk_est)');
    end

end