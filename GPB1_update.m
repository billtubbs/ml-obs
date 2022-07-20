function [xk_est,yk_est,Pk,p_seq_g_Yk] = GPB1_update(A,B,C,Q,R,T, ...
    Xkp1f_est,Ykp1f_est,Pkp1f,yk,p_seq_g_Yk)
% [xk_est,yk_est,Pk,p_seq_g_Yk] = GPB1_update(A,B,C,Q,R,T, ...
%     xkp1_est,ykp1_est,Pkp1,yk,p_seq_g_Yk)
% Update equations for simulating the generalised
% pseudo-Bayes (GPB) multi-model Kalman filter for
% state estimation of Markov jump linear systems.
%
% Based on code from the following article:
% -  Bayesian State Estimations for Markovian Jump Systems, 
%    by Shunyi Zhao, Choon Ki Ahn, Peng Shi, Yuriy S. Shmaliy, 
%    and Fei Liu, 2019.
%

    m = size(T,1);
    n = size(A{1},1);
    ny = size(C{1},1);

    updx = nan(n, m);
    updy = nan(ny, m);
    updP = nan(n, n, m);
    p_yk_g_seq_Ykm1 = nan(m, 1);
    for j = 1:m

        xkp1_est = Xkp1f_est(:,:,j);
        ykp1_est = Ykp1f_est(:,:,j);
        Pkp1 = Pkp1f(:,:,j);

        Sk = C{j} * Pkp1 * C{j}' + R{j}';
        Kf = Pkp1 * C{j}' / Sk;
        if ~isscalar(Sk)
            % Make sure covariance matrix is symmetric
            Sk = triu(Sk.',1) + tril(Sk);
        end
        p_yk_g_seq_Ykm1(j) = mvnpdf(yk, ykp1_est, Sk);
        updx(:,j) = xkp1_est + Kf * (yk - ykp1_est);
        updy(:,j) = C{j} * updx(:,j);
        updP(:,:,j) = Pkp1 - Kf * C{j} * Pkp1;
    end
    p_seq_g_Ykm1 = sum(T .* p_seq_g_Yk, 1)';
    liklihood_u = p_seq_g_Ykm1 .* p_yk_g_seq_Ykm1;
    % cf:
    % p_gamma_k = prob_gamma(gamma_k, obj.T(gamma_km1+1, :)')
    % p_seq_g_Ykm1 = p_gamma_k .* p_seq_g_Yk;
    % cond_pds = p_yk_g_seq_Ykm1 .* p_seq_g_Ykm1;
    % p_seq_g_Yk = cond_pds ./ sum(cond_pds);

    p_seq_g_Yk = liklihood_u / sum(liklihood_u);
    xk_est = sum(updx .* repmat(p_seq_g_Yk',n,1), 2);
    Pk = zeros(n,n);
    yk_est = sum(updy .* repmat(p_seq_g_Yk',ny,1), 2);
    for i = 1:m
        summP = p_seq_g_Yk(i) * (updP(:,:,i) + (updx(:,i) - xk_est) * (updx(:,i) - xk_est)');
        Pk = Pk + summP;
    end
end
