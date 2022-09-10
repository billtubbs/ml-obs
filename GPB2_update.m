function [xk_est,yk_est,Pk,Xk_est_m,Yk_est_m,Pk_m,p_seq_g_Yk] = ...
    GPB2_update(models,T,Xkp1f_est,Pkp1f,yk,p_seq_g_Yk)
% [xk_est,yk_est,Pk,Xkm_est,Ykm_est,Pkm,p_seq_g_Yk] = ...
%     GPB2_update(models,T,Xkp1f_est,Pkp1f,yk,p_seq_g_Yk)
% Update equations for simulating the second-order 
% generalised pseudo-Bayes (GPB2) multi-model Kalman 
% filter for state estimation of Markov jump linear
% systems.
%
% Steps of GPB2 algorithm at time k
%  Inputs:
%    - nj^2 prior predictions of states, x_pred_ij(k|k-1), with
%      error covariance, P_ij(k|k-1).
%    - nj^2 prior probability estimates, p(seq_ij|Y(k-1)), of the
%      nj^2 mode sequences.
%    - mode transition probabilities, T_ij.
%
%  1. Calculate the nj^2 output predictions, y_est_ij(k|k-1),
%     from the prior state predictions, x_pred_ij(k|k-1), using
%     the system models corresponding to each mode sequence.
%  2. Update the nj^2 Kalman filter state estimates, 
%     x_est_ij(k|k-1) using the current measurement, y(k), and 
%     the system models corresponding to each mode sequence.
%  3. Update the nj^2 mixing probability estimates, p(seq_ij|Y(k)), 
%     using the transition probabilities, T, the estimated 
%     likelihood of y(k) given each of the nj^2 prior output 
%     estimates, y_est_ij(k|k-1).
%  4. Use the mixing probabilities, p_seq_g_Yk to merge the
%     nj^2 state estimates into nj merged estimates, x_est_i(k|k-1),
%     for each possible system mode, i = 1, 2, ... nj, at time k.
%  5. Merge the estimates into one overall estimate, x_est(k|k-1).
%  6. Reinitialize the nj^2 Kalman filters using the nj merged
%     estimates, x_est_i(k|k-1), for the modes at time k.
%  7. Calculate the nj^2 Kalman filter predictions of the states 
%     in the next time instant, x_pred_ij(k+1|k) using the system
%     models corresponding to each mode sequence.
%

    % All models must have the same dimensions (this is not checked)
    nj = numel(models);  % number of models (modes)
    n = size(models{1}.A, 1);
    ny = size(models{1}.C, 1);
    nh = nj * nj;

    updx = nan(n, nh);
    updy = nan(ny, nh);
    updP = nan(n, n, nh);
    p_yk_g_seq_Ykm1 = nan(nh, 1);

    % Transitions modelled
    rkm1 = [1 2 1 2]';
    rk = [1 1 2 2]';

    % Transition probabilities
    p_rk_g_rkm1 = prob_transitions(rk, rkm1, T);

    for f = 1:nh

        xkp1_est = Xkp1f_est(:,:,f);
        Pkp1 = Pkp1f(:,:,f);

        % Select model according to sequence
        m = models{rk(f)};

        % Output estimate
        ykp1_est = m.C * xkp1_est;

        % Error covariance of output estimate
        Sk = m.C * Pkp1 * m.C' + m.R';

        % KF correction gain
        Kf = Pkp1 * m.C' / Sk;

        % Update Kalman filters
        updx(:,f) = xkp1_est + Kf * (yk - ykp1_est);
        updy(:,f) = m.C * updx(:,f);
        updP(:,:,f) = Pkp1 - Kf * m.C * Pkp1;

        if ~isscalar(Sk)
            % Make covariance matrix symmetric
            Sk = triu(Sk.',1) + tril(Sk);
        end
        p_yk_g_seq_Ykm1(f) = mvnpdf(yk, ykp1_est, Sk);

    end

    % Split prior probabilities from nj modes to nj^2
    p_seq_g_Yk = p_seq_g_Yk(rkm1 + 1);

    % Compute Pr(Gamma(k)|Y(k-1)) in current timestep from
    % previous estimate (Pr(Gamma(k-1)|Y(k-1))) and transition
    % probabilities
    p_seq_g_Ykm1 = p_rk_g_rkm1 .* p_seq_g_Yk;

    % Compute likelihoods of models given Y(k)
    cond_pds = p_yk_g_seq_Ykm1 .* p_seq_g_Ykm1;

    % Mixing probabilities
    p_seq_g_Yk = cond_pds / sum(cond_pds);

    % Merge the nj^2 estimates into nj modes
    [Xk_est_m, Pk_m, Yk_est_m, p_seq_g_Yk] = merge_estimates(Xkp1f_est, Pkp1f, ...
        p_seq_g_Yk, rk, nj);

    % Compute updated merged estimates
    % xk_est_out = updx * p_seq_g_Yk;
    % Above is equivalent to:
    %   sum(updx .* repmat(p_seq_g_Yk',n,1), 2)
%     yk_est_out = updy * p_seq_g_Yk;
%     Pk_out = zeros(n,n);
%     for i = 1:nj
%         Pk_out = Pk_out + p_seq_g_Yk(i) * (updP(:,:,i) + ...
%             (updx(:,i) - xk_est_out) * (updx(:,i) - xk_est_out)');
%     end

    xk_est = zeros(n, 1);
    yk_est = zeros(ny, 1);
    Pk = zeros(n, n);
    for j = 1:nh
        m_ind = rk(j);
        Xkm_est(:, m_ind) = Xkm_est(:, m_ind) + updx(:, j) .* p_seq_g_Yk(j);
        Ykm_est(:, m_ind) = Ykm_est(:, m_ind) + updy(:, j) .* p_seq_g_Yk(j);
        Pk(:,:,m_ind) = Pk(:,:,m_ind) + p_seq_g_Yk(j) * (updP(:, :, j) + ...
            (updx(:, j) - Xkm_est) * (updx(:, j) - Xkm_est)');
    end

    %[xk_est,yk_est,Pk,p_seq_g_Yk] = merge_estimates(xk_est, Pk, yk_est, p_seq_g_Yk, ...
    %r_k, nj);

end
