function xk_est = merge_estimates(xk_est, p_seq_g_Yk, gamma_k, nj)
% xk_est = merge_estimates(xk_est, p_seq_g_Yk, gamma_k, nj)
% Merges (fuses) estimates according to the modes specified
% in gamma_k.
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

end