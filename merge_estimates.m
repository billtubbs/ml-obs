function xk_est = merge_estimates(xk_est, p_seq_g_Yk, gamma_k, nj)

    n_filt = size(gamma_k, 1);
    mask = false(n_filt, nj);
    mask(sub2ind(size(mask), (1:n_filt)', gamma_k + 1)) = true;
    xk_est = sum((xk_est .* p_seq_g_Yk) .* mask, 1)';

end