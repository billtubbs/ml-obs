% Test GPB functions

clear all


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
n_filt = size(gamma_k, 1);
nj = 2;

xk_est = merge_estimates(xk_est, p_seq_g_Yk, gamma_k, nj);
assert(isequal(xk_est, [(0.3*1.1 + 0.2*1.2) (0.1*1.3 + 0.4*1.4)]'))
