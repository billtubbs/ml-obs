% Test GPB functions

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

xk_est = [1.1 1.2 1.3 1.4]';
p_seq_g_Yk = [0.3 0.2 0.1 0.4]';
gamma_k = [0 0 1 1]';
nj = 2;

xk_est = merge_estimates(xk_est, p_seq_g_Yk, gamma_k, nj);
assert(isequal(xk_est, [(0.3*1.1 + 0.2*1.2) (0.1*1.3 + 0.4*1.4)]'))


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

Pk_test = zeros(n, n);
for i = 1:nj
    summP = p_seq_g_Yk(i) * (Pkf_est(:,:,i) + (Xkf_est(:,i) - xk_est) * (Xkf_est(:,i) - xk_est)');
    Pk_test = Pk_test + summP;
end
assert(isequal(Pk, Pk_test))


%% Test mode_transitions.m

[gamma_km1, gamma_k] = mode_transitions_all(2);
assert(isequal(gamma_km1, [0 1 0 1]'))
assert(isequal(gamma_k, [0 0 1 1]'))

[gamma_km1, gamma_k] = mode_transitions_all(3);
assert(isequal(gamma_km1, [0 1 2 0 1 2 0 1 2]'))
assert(isequal(gamma_k, [0 0 0 1 1 1 2 2 2]'))
