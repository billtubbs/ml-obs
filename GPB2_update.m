function [xk_est,yk_est,Pk,p_seq_g_Yk,xk_est_out,yk_est_out,Pk_out] = ...
    GPB2_update(models,T,Xkp1f_est,Ykp1f_est,Pkp1f,yk,p_seq_g_Yk)
% [xk_est,yk_est,Pk,p_seq_g_Yk] = GPB2_update(models,T, ...
%     xkp1_est,ykp1_est,Pkp1,yk,p_seq_g_Yk)
% Update equations for simulating the second-order 
% generalised pseudo-Bayes (GPB2) multi-model Kalman 
% filter for state estimation of Markov jump linear
% systems.
%
% Steps of GPB2 algorithm at time k
%  Inputs:
%    - m^2 prior predictions of states, x_est(k|k-1),
%      error covariance, P(k|k-1), and outputs, y_est(k|k-1)
%    - prior probability estimates, p_seq_g_Ykm1, of the
%      m^2 possible mode transitions.
%    - mode transition probabilities, T
%
%  1. Update the m^2 Kalman filters using the current 
%     measurement, yk, and the current models.
%  2. Update the m^2 mixing probability estimates, p_seq_g_Yk, 
%     using the m^2 transition probabilities, T, the estimated 
%     likelihood of yk given each of the m^2 prior output 
%     estimates, y_est(k|k-1).
%  3. Use the mixing probabilities, p_seq_g_Yk to merge the
%     m^2 estimates into m estimates for each possible mode
%     at time k.
%  4. Merge the estimates into one overall estimate at time
%     k.
%  5. Reinitialize the m^2 Kalman filters using the m
%     estimates for the possible modes at time k.
%  6. Run the m^2 Kalman filter predictions using the 
%     the current models to calculate prior estimates of the
%     states at time k + 1.
%

    nj = numel(models);  % number of models (modes)
    % All models must be of same dimensions
    n = size(models{1}.A, 1);
    ny = size(models{1}.C, 1);
    n_filt = nj * nj;

    updx = nan(n, n_filt);
    updy = nan(ny, n_filt);
    updP = nan(n, n, n_filt);
    p_yk_g_seq_Ykm1 = nan(n_filt, 1);

    % Transitions modelled
    gamma_km1 = [0 1 0 1]';
    gamma_k = [0 0 1 1]';

    % Transition probabilities
    p_gamma_k_g_gamma_km1 = prob_transitions(gamma_k, gamma_km1, T);

    for f = 1:n_filt

        xkp1_est = Xkp1f_est(:,:,f);
        ykp1_est = Ykp1f_est(:,:,f);
        Pkp1 = Pkp1f(:,:,f);

        % Select model according to sequence
        m = models{gamma_k(f) + 1};  % TODO: Change to gamma_km1 later

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

    % Split prior probabilities from m modes to m^2
    p_seq_g_Yk = p_seq_g_Yk(gamma_km1 + 1);

    % Compute Pr(Gamma(k)|Y(k-1)) in current timestep from
    % previous estimate (Pr(Gamma(k-1)|Y(k-1))) and transition
    % probabilities
    p_seq_g_Ykm1 = p_gamma_k_g_gamma_km1 .* p_seq_g_Yk;

    % Compute likelihoods of models given Y(k)
    cond_pds = p_yk_g_seq_Ykm1 .* p_seq_g_Ykm1;

    % Mixing probabilities
    mask = false(n_filt, nj);
    mask(sub2ind(size(mask), (1:n_filt)', gamma_k + 1)) = true;

    p_mix = Mix_u(:,j) / sum(Mix_u(:,j))
    p_seq_g_Yk = cond_pds / sum(cond_pds);

    % Create masks to merge the m^2 estimates into m modes
    %mask = false(n_filt, nj);
    %mask(sub2ind(size(mask), (1:n_filt)', gamma_k + 1)) = true;

    % Compute updated merged estimates
    xk_est_out = updx * p_seq_g_Yk;
    % Above is equivalent to:
    %   sum(updx .* repmat(p_seq_g_Yk',n,1), 2)
    yk_est_out = updy * p_seq_g_Yk;
    Pk_out = zeros(n,n);
    for i = 1:nj
        Pk_out = Pk_out + p_seq_g_Yk(i) * (updP(:,:,i) + ...
            (updx(:,i) - xk_est_out) * (updx(:,i) - xk_est_out)');
    end

    % Compute updated semi-merged estimates
    %xk_est = sum((updx' .* p_seq_g_Yk) .* mask, 1);
    %yk_est = sum((updy' .* p_seq_g_Yk) .* mask, 1);
    xk_est = zeros(n, nj);
    yk_est = zeros(ny, nj);
    Pk = zeros(n, n, nj);
    %xk_est_out = zeros(n, 1);
    %yk_est_out = zeros(ny, 1);
    %Pk_out = zeros(n, n);
    for j = 1:n_filt
        m_ind = gamma_k(j) + 1;
        xk_est(:, m_ind) = xk_est(:, m_ind) + updx(:, j) .* p_seq_g_Yk(j);
        yk_est(:, m_ind) = yk_est(:, m_ind) + updy(:, j) .* p_seq_g_Yk(j);
        Pk(:,:,m_ind) = Pk(:,:,m_ind) + p_seq_g_Yk(j) * (updP(:, :, j) + ...
            (updx(:, j) - xk_est) * (updx(:, j) - xk_est)');
    end

    %[xk_est,yk_est,Pk,p_seq_g_Yk] = merge_estimates(xk_est, Pk, yk_est, p_seq_g_Yk, ...
    %gamma_k, nj);

end
