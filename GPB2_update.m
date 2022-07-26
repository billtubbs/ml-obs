function [xk_est,yk_est,Pk,p_seq_g_Yk,xk_est_out,yk_est_out,Pk_out] = ...
    GPB2_update(A,B,C,Q,R,T,Xkp1f_est,Ykp1f_est,Pkp1f,yk,p_seq_g_Yk)
% [xk_est,yk_est,Pk,p_seq_g_Yk] = GPB2_update(A,B,C,Q,R,T, ...
%     xkp1_est,ykp1_est,Pkp1,yk,p_seq_g_Yk)
% Update equations for simulating the second-order 
% generalised pseudo-Bayes (GPB2) multi-model Kalman 
% filter for state estimation of Markov jump linear
% systems.
%
% Based on code from the following article:
% -  Bayesian State Estimations for Markovian Jump Systems, 
%    by Shunyi Zhao, Choon Ki Ahn, Peng Shi, Yuriy S. Shmaliy, 
%    and Fei Liu, 2019.
%

    nj = numel(A);  % number of models (modes)
    n_filt = nj * nj;  % number of hypotheses (filters)
    n = size(A{1},1);  % number of states
    ny = size(C{1},1);  % number of outputs

    updx = nan(n, n_filt);
    updy = nan(ny, n_filt);
    updP = nan(n, n, n_filt);
    p_yk_g_seq_Ykm1 = nan(n_filt, 1);

    % Transitions modelled
    gamma_km1 = [0 1 0 1]';
    gamma_k = [0 0 1 1]';

    % Transition probabilities
    p_gamma_k_g_gamma_km1 = prob_transitions(gamma_k, gamma_km1, T);

    for j = 1:n_filt

        xkp1_est = Xkp1f_est(:,:,j);
        ykp1_est = Ykp1f_est(:,:,j);
        Pkp1 = Pkp1f(:,:,j);

        % Select model according to sequence
        m_ind = gamma_k(j) + 1;
        Sk = C{m_ind} * Pkp1 * C{m_ind}' + R{m_ind}';
        Kf = Pkp1 * C{m_ind}' / Sk;
        if ~isscalar(Sk)
            % Make sure covariance matrix is symmetric
            Sk = triu(Sk.',1) + tril(Sk);
        end
        p_yk_g_seq_Ykm1(j) = mvnpdf(yk, ykp1_est, Sk);
        updx(:,j) = xkp1_est + Kf * (yk - ykp1_est);
        updy(:,j) = C{m_ind} * updx(:,j);
        updP(:,:,j) = Pkp1 - Kf * C{m_ind} * Pkp1;

    end

    % Split prior probabilities from m modes to m^2
    p_seq_g_Yk = p_seq_g_Yk(gamma_k + 1);

    % Compute Pr(Gamma(k)|Y(k-1)) in current timestep from
    % previous estimate (Pr(Gamma(k-1)|Y(k-1))) and transition
    % probabilities
    p_seq_g_Ykm1 = p_gamma_k_g_gamma_km1 .* p_seq_g_Yk;

    % Compute likelihoods of models given Y(k)
    cond_pds = p_seq_g_Ykm1 .* p_yk_g_seq_Ykm1;
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
            (updx(:,i) - xk_est) * (updx(:,i) - xk_est)');
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

    p_seq_g_Yk = merge_estimates(p_seq_g_Yk, p_seq_g_Yk, gamma_k, nj);

end
