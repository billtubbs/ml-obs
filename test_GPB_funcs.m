% Test the following GPB functions:
%  - prob_transitions.m
%  - merge_estimates.m
%  - weighted_avg_estimates.m
%  - mode_transitions_all.m
%

clear all
rng(0)


%% Test prob_transitions.m

% Transition probability matrix (2 modes)
T = [0.95 0.05; 0.01 0.99];

% Transitions (x4)
gamma_km1 = [0 1 0 1]';
gamma_k = [0 0 1 1]';
n_filt = size(gamma_k, 1);

% Calculate transition probabilities
p_gamma_k_g_gamma_km1 = prob_transitions(gamma_k, gamma_km1, T);
assert(isequal(p_gamma_k_g_gamma_km1, reshape(T, [], 1)));

% Simple way to calculate them
p_test = zeros(n_filt, 1);
for i = 1:n_filt
    p_test(i) = T(gamma_km1(i) + 1, gamma_k(i) + 1);
end
assert(isequal(p_gamma_k_g_gamma_km1, p_test))

% Alternative way to calculate them
idx = sub2ind(size(T), gamma_km1+1, gamma_k+1);
assert(isequal(idx, (1:4)'))
assert(isequal(p_gamma_k_g_gamma_km1, T(idx)))


%% Test merge_estimates.m

% 1. Single state and output
xk_est = [1.1 1.2 1.3 1.4];
Pk = cat(3,100,110,120,130);
yk_est = [2 2.2 2.6 2.8];
p_seq_g_Yk = [0.3 0.2 0.1 0.4]';
gamma_k = [0 0 1 1]';
nj = 2;

% Function to test
[xk_est_m,yk_est_m,Pk_m,p_seq_g_Yk_m] = ...
    merge_estimates(xk_est, Pk, yk_est, p_seq_g_Yk, gamma_k, nj);

% Test outputs
xk_est_m_test = [(0.3*1.1+0.2*1.2)/0.5 (0.1*1.3+0.4*1.4)/0.5];
assert(isequal(xk_est_m, xk_est_m_test))
yk_est_m_test = [(0.3*2+0.2*2.2)/0.5 (0.1*2.6+0.4*2.8)/0.5];
assert(isequal(yk_est_m, yk_est_m_test))
Pk_m_test = cat(3, 0, 0);
Pk_m_test(:,:,1) = ...
      p_seq_g_Yk(1) * (Pk(:,:,1) + (1.1 - xk_est_m_test(1)).^2) / 0.5 ...
    + p_seq_g_Yk(2) * (Pk(:,:,2) + (1.2 - xk_est_m_test(1)).^2) / 0.5;
Pk_m_test(:,:,2) = ...
      p_seq_g_Yk(3) * (Pk(:,:,3) + (1.3 - xk_est_m_test(2)).^2) / 0.5 ...
    + p_seq_g_Yk(4) * (Pk(:,:,4) + (1.4 - xk_est_m_test(2)).^2) / 0.5;
assert(isequal(Pk_m, Pk_m_test))
assert(isequal(p_seq_g_Yk_m, [0.5; 0.5]))

% 2. Multiple states and outputs
xk_est = [1.1 1.2 1.3 1.4; ...
          2.1 2.2 2.3 2.4; ...
          3.1 3.2 3.3 3.4];
n = size(xk_est, 1);
Pk = cat(3,10*eye(n),20*eye(n),30*eye(n),40*eye(n));
yk_est = [0.5 1 0;
          0.5 0 1] * xk_est;
p_seq_g_Yk = [0.3 0.2 0.1 0.4]';
gamma_k = [0 0 1 1]';
nj = 2;

% Function to test
[xk_est_m,yk_est_m,Pk_m,p_seq_g_Yk_m] = ...
    merge_estimates(xk_est, Pk, yk_est, p_seq_g_Yk, gamma_k, nj);

% Compare with this code adapted from GPB2_estimation.m
Mix_u = reshape(p_seq_g_Yk, nj, nj);
for j = 1:nj

    % Normalize the mode probabilities
    Updx = xk_est(:,(gamma_k + 1) == j);
    Upd_u = Mix_u(:,j) / sum(Mix_u(:,j));
    UpdP = Pk(:,:,(gamma_k + 1) == j);

    % Mix the estimates
    Mix_x(:,j) = sum(Updx .* repmat(Upd_u', n, 1), 2);
    Mix_P(:,:,j) = zeros(n, n);
    for i = 1:nj
        summP = Upd_u(i) * (UpdP(:,:,i) + ...
            (Updx(:,i) - Mix_x(:,j)) * (Updx(:,i) - Mix_x(:,j))');
        Mix_P(:,:,j) =  Mix_P(:,:,j) + summP;
    end

    Storu(j) = sum(Mix_u(:,j));

end

% Tests
assert(isequal(Mix_x(:,1), ...
    sum(xk_est(:, 1:2) .* p_seq_g_Yk(1:2)', 2) / sum(p_seq_g_Yk(1:2))))
assert(isequal(Mix_x(:,2), ...
    sum(xk_est(:, 3:4) .* p_seq_g_Yk(3:4)', 2) / sum(p_seq_g_Yk(3:4))))
% Above is equivalent to
assert(isequal(Mix_x(:,1), ...
    sum(xk_est(:, gamma_k == 0) .* p_seq_g_Yk(gamma_k == 0)', 2) / sum(p_seq_g_Yk(gamma_k == 0))))
assert(isequal(Mix_x(:,2), ...
    sum(xk_est(:, gamma_k == 1) .* p_seq_g_Yk(gamma_k == 1)', 2) / sum(p_seq_g_Yk(gamma_k == 1))))

assert(isequal(xk_est_m, Mix_x));
assert(isequal(Pk_m, Mix_P))


%% Test weighted_avg_estimates.m

nj = 3;
n = 4;
ny = 2;
Xkf_est = randn(n, 1, nj);
Ykf_est = randn(ny, 1, nj);
Pkf_est = randn(n, n, nj);

p_seq_g_Yk = randn(nj, 1);
p_seq_g_Yk = p_seq_g_Yk ./ sum(p_seq_g_Yk);

[xk_est, yk_est, Pk] = weighted_avg_estimates(Xkf_est, Ykf_est, ...
    Pkf_est, p_seq_g_Yk);

assert(isequal(xk_est, Xkf_est(:,:,1) * p_seq_g_Yk(1) + ...
    Xkf_est(:,:,2) * p_seq_g_Yk(2) + Xkf_est(:,:,3) * p_seq_g_Yk(3)))
assert(isequal(yk_est, Ykf_est(:,:,1) * p_seq_g_Yk(1) + ...
    Ykf_est(:,:,2) * p_seq_g_Yk(2) + Ykf_est(:,:,3) * p_seq_g_Yk(3)))

Pk_m_test = zeros(n, n);
for i = 1:nj
    summP = p_seq_g_Yk(i) * (Pkf_est(:,:,i) + (Xkf_est(:,i) - xk_est) * (Xkf_est(:,i) - xk_est)');
    Pk_m_test = Pk_m_test + summP;
end
assert(isequal(Pk, Pk_m_test))


%% Test mode_transitions_all.m

[gamma_km1, gamma_k] = mode_transitions_all(2);
assert(isequal(gamma_km1, [0 1 0 1]'))
assert(isequal(gamma_k, [0 0 1 1]'))

[gamma_km1, gamma_k] = mode_transitions_all(3);
assert(isequal(gamma_km1, [0 1 2 0 1 2 0 1 2]'))
assert(isequal(gamma_k, [0 0 0 1 1 1 2 2 2]'))
