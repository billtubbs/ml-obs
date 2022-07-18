% Update equations for simulating the generalised
% pseudo-Bayes (GPB) multi-model Kalman filter for
% state estimation of Markov jump linear systems.
%
% Update equations for generalised pseudo-Bayes multi-
% model observer algorithm GBP1.
%
% Based on code from the following article:
% -  Bayesian State Estimations for Markovian Jump Systems, 
%    by Shunyi Zhao, Choon Ki Ahn, Peng Shi, Yuriy S. Shmaliy, 
%    and Fei Liu, 2019.
%

function [xk_est,yk_est,Pk,p_seq_g_Yk] = GPB1_update(A,B,C,Q,R,T,xk_est,Pk,yk,uk,p_seq_g_Yk)
    m = size(T,1);
    n = size(A{1},1);
    ny = size(C{1},1);
    liklihood_u = nan(m, 1);
    updx = nan(n, m);
    updy = nan(ny, m);
    updP = nan(n, n, m);
    for j = 1:m
        Pred_u = T(:,j) .* p_seq_g_Yk;
        xkp1_est = A{j} * xk_est + B{j} * uk;
        Pkp1 = A{j} * Pk * A{j}' + Q{j};
        ykp1_est = C{j} * xkp1_est;
        S = C{j} * Pkp1 * C{j}' + R{j}';
        K = Pkp1 * C{j}' / S;
        if ~isscalar(S)
            % Make sure covariance matrix is symmetric
            S = triu(S.',1) + tril(S);
        end
        liklihood_u(j) = sum(Pred_u) * mvnpdf(yk, ykp1_est, S);
        updx(:,j) = xkp1_est + K * (yk - ykp1_est);
        updy(:,j) = C{j} * updx(:,j);
        updP(:,:,j) = Pkp1 - K * C{j} * Pkp1;
    end
    p_seq_g_Yk = liklihood_u / sum(liklihood_u);
    xk_est = sum(updx .* repmat(p_seq_g_Yk',n,1), 2);
    Pk = zeros(n,n);
    yk_est = sum(updy .* repmat(p_seq_g_Yk',ny,1), 2);
    for i = 1:m
        summP = p_seq_g_Yk(i) * (updP(:,:,i) + (updx(:,i) - xk_est) * (updx(:,i) - xk_est)');
        Pk = Pk + summP;
    end
end
