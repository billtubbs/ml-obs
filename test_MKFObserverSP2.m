% Test functions mkf_observer_AFMM.m and update_AFMM.m

clear all
%plot_dir = 'plots';

seed = 0;
rng(seed)


%% Test initialization with rodin_step system

% Load system and disturbance model from file
sys_rodin_step

% Load observers from file
obs_rodin_step

% Observer model without disturbance noise input
Bu = B(:, u_meas);
Du = D(:, u_meas);
nu = sum(u_meas);
nw = sum(~u_meas);

% Set noise variances for observer design
sigma_M = 0.1;
sigma_W = [0; 0];

% Check observer attributes
assert(strcmp(MKF_SP21.type, "MKF-SP"))
assert(MKF_SP21.epsilon == 0.01)
assert(isequal(MKF_SP21.sigma_wp, sigma_wp))
assert(MKF_SP21.n_filt == 5)
assert(MKF_SP21.n_min == 3)
assert(isequal(MKF_SP21.n_hold, 3))
assert(isequal(MKF_SP21.n_main, 2))
assert(isequaln(MKF_SP21.f_hold, [3 4 5]))
assert(isequaln(MKF_SP21.f_main, [1 2]))
assert(isequaln(MKF_SP21.i, [0 0]))
assert(MKF_SP21.n == 2)
assert(MKF_SP21.nu == 1)
assert(MKF_SP21.ny == 1)
assert(MKF_SP21.nj == 2)
assert(isequal(MKF_SP21.A{1}, A) && isequal(MKF_SP21.A{2}, A))
assert(isequal(MKF_SP21.B{1}, Bu) && isequal(MKF_SP21.B{2}, Bu))
assert(isequal(MKF_SP21.C{1}, C) && isequal(MKF_SP21.C{2}, C))
assert(isequal(MKF_SP21.D{1}, Du) && isequal(MKF_SP21.D{2}, Du))
assert(MKF_SP21.Ts == Ts)
assert(isequaln(MKF_SP21.u_meas, u_meas))
assert(isequal(MKF_SP21.Q{1}, [0.01 0; 0 sigma_wp(1)^2]))
assert(isequal(MKF_SP21.Q{2}, [0.01 0; 0 sigma_wp(2)^2]))
assert(isequal(MKF_SP21.R{1}, R) && isequal(MKF_SP21.R{2}, R))
assert(numel(MKF_SP21.filters) == MKF_SP21.n_filt)
assert(isequal(size(MKF_SP21.seq), [MKF_SP21.n_filt 1]))
assert(isequal(size(cell2mat(MKF_SP21.seq)), [MKF_SP21.n_filt MKF_SP21.f]))
assert(MKF_SP21.f == size(MKF_SP21.seq{1}, 2))
assert(isequal(size(MKF_SP21.xkp1_est), [n 1]))
assert(isequal(size(MKF_SP21.ykp1_est), [ny 1]))
assert(isequal(MKF_SP21.p_gamma, [1-MKF_SP21.epsilon; MKF_SP21.epsilon]))

% Check initialization of filters
assert(isequal(MKF_SP21.filters{1}.P, MKF_SP21.P0))
assert(isequal(MKF_SP21.filters{1}.x0, MKF_SP21.x0))
for i = 2:MKF_SP21.n_filt
    assert(isequal(MKF_SP21.filters{i}.P, 1e10*eye(2)))
    assert(isequal(MKF_SP21.filters{i}.x0, MKF_SP21.x0))
end

assert(MKF_SP22.epsilon == 0.01)
assert(isequal(MKF_SP22.sigma_wp, sigma_wp))
assert(MKF_SP22.n_filt == 10)
assert(MKF_SP22.n_min == 4)
assert(isequal(MKF_SP22.n_hold, 4))
assert(isequal(MKF_SP22.n_main, 6))
assert(isequaln(MKF_SP22.f_hold, [7 8 9 10]))
assert(isequaln(MKF_SP22.f_main, [1 2 3 4 5 6]))
assert(isequaln(MKF_SP22.i, [0 0]))
assert(MKF_SP22.n == 2)
assert(MKF_SP22.nu == 1)
assert(MKF_SP22.ny == 1)
assert(MKF_SP22.nj == 2)
assert(isequal(MKF_SP22.A{1}, A) && isequal(MKF_SP22.A{2}, A))
assert(isequal(MKF_SP22.B{1}, Bu) && isequal(MKF_SP22.B{2}, Bu))
assert(isequal(MKF_SP22.C{1}, C) && isequal(MKF_SP22.C{2}, C))
assert(isequal(MKF_SP22.D{1}, Du) && isequal(MKF_SP22.D{2}, Du))
assert(isequal(MKF_SP22.B{1}, Bu) && isequal(MKF_SP22.B{2}, Bu))
assert(isequal(MKF_SP22.C{1}, C) && isequal(MKF_SP22.C{2}, C))
assert(isequal(MKF_SP22.D{1}, Du) && isequal(MKF_SP22.D{2}, Du))
assert(MKF_SP22.Ts == Ts)
assert(isequaln(MKF_SP22.u_meas, u_meas))
assert(isequal(MKF_SP22.Q{1}, [0.01 0; 0 sigma_wp(1)^2]))
assert(isequal(MKF_SP22.Q{2}, [0.01 0; 0 sigma_wp(2)^2]))
assert(isequal(MKF_SP22.R{1}, R) && isequal(MKF_SP22.R{2}, R))
assert(numel(MKF_SP22.filters) == MKF_SP22.n_filt)
assert(isequal(size(MKF_SP22.seq), [MKF_SP22.n_filt 1]))
assert(isequal(size(cell2mat(MKF_SP22.seq)), [MKF_SP22.n_filt MKF_SP22.f]))
assert(MKF_SP22.f == size(MKF_SP22.seq{1}, 2))
assert(isequal(size(MKF_SP22.xkp1_est), [n 1]))
assert(isequal(size(MKF_SP22.ykp1_est), [ny 1]))
assert(isequal(MKF_SP22.p_gamma, [1-MKF_SP22.epsilon; MKF_SP22.epsilon]))

% Check initialization of filters
assert(isequal(MKF_SP22.filters{1}.P, MKF_SP22.P0))
assert(isequal(MKF_SP22.filters{1}.x0, MKF_SP22.x0))
for i = 2:MKF_SP22.n_filt
    assert(isequal(MKF_SP22.filters{i}.P, 1e10*eye(2)))
    assert(isequal(MKF_SP22.filters{i}.x0, MKF_SP22.x0))
end

% Check optional definition with an initial state estimate works
x0 = [0.1; 0.5];
MKF_SP_testx0 = MKFObserverSP2(A,B,C,D,Ts,u_meas,P0,epsilon,sigma_wp, ...
    Q0,R,n_filt,f,n_min,label,x0);
assert(isequal(MKF_SP_testx0.xkp1_est, x0))
assert(isequal(MKF_SP_testx0.ykp1_est, C * x0))
for i = 1:MKF_SP_testx0.n_filt
    assert(isequal(MKF_SP_testx0.filters{i}.x0, MKF_SP_testx0.x0))
end


%% Test convergence to steady-state

% Load system and disturbance model from file
sys_rodin_step

% Load observers from file
obs_rodin_step

% Check steady-state at x0 = [0; 0]
obs = MKF_SP21;
assert(isequal(obs.xkp1_est, [0; 0]))
assert(isequal(obs.ykp1_est, 0))
nT = 10;
U_m = zeros(nT+1, sum(u_meas));
Y_m = zeros(nT+1, ny);
for i = 1:(nT+1)
    uk = U_m(i,:)';
    yk = Y_m(i,:)';
    obs.update(yk, uk);
    assert(isequal(obs.xkp1_est, [0; 0]))
    assert(isequal(obs.ykp1_est, 0))
end

% Check steady-state at x0 = [1; 0]
x0 = [1; 0];
obs = MKFObserverSP2(A,B,C,D,Ts,u_meas,P0,epsilon,sigma_wp, ...
    Q0,R,n_filt,f,n_min,label,x0);
assert(isequal(obs.xkp1_est, x0))
assert(isequal(obs.ykp1_est, 0.3))
nT = 10;
U_m = 0.3*ones(nT+1, 1);
U = [U_m zeros(nT+1,1)];
t = Ts*(0:nT)';
[Y,~,X] = lsim(Gpss,U,t,x0);
assert(all(Y == Y(1,1)))
for i = 1:(nT+1)
    uk = U_m(i,:)';
    yk = Y(i,:)';
    obs.update(yk, uk);
    assert(all(abs(obs.xkp1_est - x0) < 1e-6))
    assert(abs(obs.ykp1_est - 0.3) < 1e-6)
end


%% Test sequence updates on SISO linear system

% Load system and disturbance model from file
sys_rodin_step

% Load observers from file
obs_rodin_step

nT = 6;
x0 = [0; 0];
n_filt = 5;
f = 8;
n_min = 2;
obs = MKFObserverSP2(A,B,C,D,Ts,u_meas,P0,epsilon,sigma_wp, ...
    Q0,R,n_filt,f,n_min,label,x0);
assert(isequal(obs.xkp1_est, x0))
assert(isequal(obs.ykp1_est, 0))
assert(obs.d == 1)

% % Generate test simulation data
% nT = 10;
% U_m = zeros(nT+1, sum(u_meas));
% % Add a random shock
% Wp = zeros(nT+1, sum(~u_meas));
% Wp(5, :) = 1;
% % Compute outputs (n0 measurement noise)
% [Y_m, t] = lsim(Gpss, [U_m Wp], t);
% [t U_m Wp Y_m]

% Test simuation data
sim_data = [ ...
         0         0         0         0;
    0.5000         0         0         0;
    1.0000         0         0         0;
    1.5000         0         0         0;
    2.0000         0    1.0000         0;
    2.5000         0         0         0;
    3.0000         0         0    0.3000;
    3.5000         0         0    0.5100;
    4.0000         0         0    0.6570;
    4.5000         0         0    0.7599;
    5.0000         0         0    0.8319];
nT = size(sim_data, 1);
t = sim_data(:, 1);
U_m = sim_data(:, 2);
Wp = sim_data(:, 3);
Y_m = sim_data(:, 4);

% Set marker values on each sequence - for testing only
% these values at the end of the sequences are not used
% by the observer.
for i = 1:n_filt
    obs.seq{i}(8) = i;
end
seq0 = [
    0 0 0 0 0 0 0 1
    0 0 0 0 0 0 0 2
    0 0 0 0 0 0 0 3
    0 0 0 0 0 0 0 4
    0 0 0 0 0 0 0 5
];
%%disp(obs.i)  % use for debugging
%disp(debug_array(obs))
assert(isequaln(obs.i, [0 0]))
assert(isequaln(obs.i_next, [1 1]))
assert(isequaln(cell2mat(obs.seq), seq0))
assert(isequal(obs.n_hold, 2))
assert(isequal(obs.n_main, 3))
assert(isequaln(obs.f_hold, [4 5]))
assert(isequaln(obs.f_main, [1 2 3]))

% Check probabilities
assert(isequal(obs.gamma_k, [0 0 0 0 0]'))
assert(isequal(obs.p_gamma_k, [0 0 0 0 0]'))
assert(isequal(obs.p_yk_g_seq_Ykm1, [0 0 0 0 0]'))
assert(isequal(obs.p_seq_g_Ykm1, [0 0 0 0 0]'))
assert(isequal(obs.p_seq_g_Yk, [1 zeros(1, 4)]'))

% Check initialization of filters
assert(isequal(obs.filters{1}.P, obs.P0))
assert(isequal(obs.filters{2}.P, 1e10*eye(2)))
assert(isequal(obs.filters{3}.P, 1e10*eye(2)))
assert(isequal(obs.filters{4}.P, 1e10*eye(2)))
assert(isequal(obs.filters{5}.P, 1e10*eye(2)))

% Check estimates
assert(isequal(obs.xkp1_est, [0; 0]))
assert(isequal(obs.ykp1_est, 0))

% Update at k = 0
i = 1;
uk = U_m(i,:)';
yk = Y_m(i,:)';
obs.update(yk, uk);
seq = [
    0 0 0 0 0 0 0 1
    1 0 0 0 0 0 0 1
    0 0 0 0 0 0 0 3
    0 0 0 0 0 0 0 4
    0 0 0 0 0 0 0 5
];
%disp(obs.i)
%disp(debug_array(obs))
assert(isequaln(obs.i, [1 1]))
assert(isequaln(obs.i_next, [2 1]))
assert(isequaln(cell2mat(obs.seq), seq))
assert(isequaln(obs.f_hold, [2 4]))
assert(isequaln(obs.f_main, [1 5 3]))

% Check probabilities
assert(isequal(obs.gamma_k, [0 1 0 0 0]'))
assert(isequal(obs.p_gamma_k, [0.99 0.01 0.99 0.99 0.99]'))
assert(isequal(round(obs.p_yk_g_seq_Ykm1, 4), [0.0420 0.0420 0 0 0]'))
assert(isequal(round(obs.p_seq_g_Ykm1, 4), [0.99 0.01 0 0 0]'))
assert(isequal(round(obs.p_seq_g_Yk, 4), [0.99 0.01 0 0 0]'))

% Check estimates
assert(isequal(obs.xkp1_est, [0; 0]))
assert(isequal(obs.ykp1_est, 0))
assert(isequal(obs.ykp1_est - yk, 0))

% Update at k = 1
i = 2;
uk = U_m(i,:)';
yk = Y_m(i,:)';
obs.update(yk, uk);
seq = [
    0 0 0 0 0 0 0 1
    1 0 0 0 0 0 0 1
    0 0 0 0 0 0 0 3
    0 0 0 0 0 0 0 4
    0 1 0 0 0 0 0 1
];
%disp(obs.i)
%disp(debug_array(obs))
assert(isequaln(obs.i, [2 1]))
assert(isequaln(obs.i_next, [3 1]))
assert(isequaln(cell2mat(obs.seq), seq))
assert(isequaln(obs.f_hold, [5 2]))
assert(isequaln(obs.f_main, [1 4 3]))

% Check probabilities
assert(isequal(obs.gamma_k, [0 0 0 0 1]'))
assert(isequal(obs.p_gamma_k, [0.99 0.99 0.99 0.99 0.01]'))
assert(isequal(round(obs.p_yk_g_seq_Ykm1, 4), [0.0420 0.0420 0 0 0.0420]'))
assert(isequal(round(obs.p_seq_g_Ykm1, 4), [0.9801 0.0099 0 0 0.0099]'))  % NOTE: doesn't quite sum to 1
assert(isequal(round(obs.p_seq_g_Yk, 4), [0.9802 0.0099 0 0 0.0099]'))

% Check estimates
assert(isequal(obs.xkp1_est, [0; 0]))
assert(isequal(obs.ykp1_est, 0))
assert(isequal(obs.ykp1_est - yk, 0))

% Update at k = 2
i = 3;
uk = U_m(i,:)';
yk = Y_m(i,:)';
obs.update(yk, uk);
seq = [
    0 0 0 0 0 0 0 1
    1 0 0 0 0 0 0 1
    0 0 0 0 0 0 0 3
    0 0 1 0 0 0 0 1
    0 1 0 0 0 0 0 1
];
%disp(obs.i)
%disp(debug_array(obs))
assert(isequaln(obs.i, [3 1]))
assert(isequaln(obs.i_next, [4 1]))
assert(isequaln(cell2mat(obs.seq), seq))
assert(isequaln(obs.f_hold, [4 5]))
assert(isequaln(obs.f_main, [1 2 3]))

% Check probabilities
assert(isequal(obs.gamma_k, [0 0 0 1 0]'))
assert(isequal(obs.p_gamma_k, [0.99 0.99 0.99 0.01 0.99]'))
assert(isequal(round(obs.p_yk_g_seq_Ykm1, 4), [1.8682 1.0834 1.8680 1.8682 1.8682]'))
assert(isequal(round(obs.p_seq_g_Ykm1, 4), [0.9704 0.0098 0 0.0098 0.0098]'))  % NOTE: doesn't quite sum to 1
assert(isequal(round(obs.p_seq_g_Yk, 4), [0.9746 0.0057 0 0.0098 0.0098]'))

% Check estimates
assert(isequal(obs.xkp1_est, [0; 0]))
assert(isequal(obs.ykp1_est, 0))
assert(isequal(obs.ykp1_est - yk, 0))

% Update at k = 3
i = 4;
uk = U_m(i,:)';
yk = Y_m(i,:)';
obs.update(yk, uk);
seq = [
    0 0 0 0 0 0 0 1
    1 0 0 0 0 0 0 1
    0 0 0 1 0 0 0 1
    0 0 1 0 0 0 0 1
    0 1 0 0 0 0 0 1  % all sequences are now splits from #1
];
%disp(obs.i)
%disp(debug_array(obs))
assert(isequaln(obs.i, [4 1]))
assert(isequaln(obs.i_next, [5 1]))
assert(isequaln(cell2mat(obs.seq), seq))
assert(isequaln(obs.f_hold, [3 4]))
assert(isequaln(obs.f_main, [1 2 5]))

% Check probabilities
assert(isequal(obs.gamma_k, [0 0 1 0 0]'))
assert(isequal(obs.p_gamma_k, [0.99 0.99 0.01 0.99 0.99]'))
assert(isequal(round(obs.p_yk_g_seq_Ykm1, 4), [2.4676 2.0186 2.4676 2.4676 1.1707]'))
assert(isequal(round(obs.p_seq_g_Ykm1, 4), [0.9649 0.0057 0.0097 0.0097 0.0097]'))
assert(isequal(round(obs.p_seq_g_Yk, 4), [0.9711 0.0047 0.0098 0.0098 0.0047]'))

% Check estimates
assert(isequal(obs.xkp1_est, [0; 0]))
assert(isequal(obs.ykp1_est, 0))
assert(isequal(obs.ykp1_est - yk, 0))

% Update at k = 4
i = 5;
uk = U_m(i,:)';
yk = Y_m(i,:)';
obs.update(yk, uk);
seq = [
    0 0 0 0 0 0 0 1
    0 0 0 0 1 0 0 1  % this seq. has now been replaced
    0 0 0 1 0 0 0 1
    0 0 1 0 0 0 0 1
    0 1 0 0 0 0 0 1 
];
%disp(obs.i)
%disp(debug_array(obs))
assert(isequaln(obs.i, [5 1]))
assert(isequaln(obs.i_next, [6 1]))
assert(isequaln(cell2mat(obs.seq), seq))
assert(isequaln(obs.f_hold, [2 3]))
assert(isequaln(obs.f_main, [1 4 5]))

% Check probabilities
assert(isequal(obs.gamma_k, [0 1 0 0 0]'))
assert(isequal(obs.p_gamma_k, [0.99 0.01 0.99 0.99 0.99]'))
assert(isequal(round(obs.p_yk_g_seq_Ykm1, 4), [2.8078 2.8078 2.8078 1.2019 2.0187]'))
assert(isequal(round(obs.p_seq_g_Ykm1, 4), [0.9614 0.0097 0.0097 0.0097 0.0046]'))
assert(isequal(round(obs.p_seq_g_Yk, 4), [0.9728 0.0098 0.0098 0.0042 0.0034]'))

% Check estimates
assert(isequal(obs.xkp1_est, [0; 0]))
assert(isequal(obs.ykp1_est, 0))
assert(isequal(obs.ykp1_est - yk, 0))

% Update at k = 5
i = 6;
uk = U_m(i,:)';
yk = Y_m(i,:)';
obs.update(yk, uk);
seq = [
    0 0 0 0 0 0 0 1
    0 0 0 0 1 0 0 1
    0 0 0 1 0 0 0 1
    0 0 1 0 0 0 0 1
    0 0 0 0 0 1 0 1
];
%disp(obs.i)
%disp(debug_array(obs))
assert(isequaln(obs.i, [6 1]))
assert(isequaln(obs.i_next, [7 1]))
assert(isequaln(cell2mat(obs.seq), seq))
assert(isequaln(obs.f_hold, [5 2]))
assert(isequaln(obs.f_main, [1 4 3]))

% Check probabilities
assert(isequal(obs.gamma_k, [0 0 0 0 1]'))
assert(isequal(obs.p_gamma_k, [0.99 0.99 0.99 0.99 0.01]'))
assert(isequal(round(obs.p_yk_g_seq_Ykm1, 4), [3.0218 3.0218 1.2172 2.0223 3.0218]'))
assert(isequal(round(obs.p_seq_g_Ykm1, 4), [0.9631 0.0097 0.0097 0.0042 0.0097]'))
assert(isequal(round(obs.p_seq_g_Yk, 4), [0.9736 0.0098 0.0040 0.0028 0.0098]'))

% Check estimates
assert(isequal(obs.xkp1_est, [0; 0]))
assert(isequal(obs.ykp1_est, 0))
assert(isequal(obs.ykp1_est - yk, 0))

% Update at k = 6  *** First non-zero measurement ***
i = 7;
uk = U_m(i,:)';
yk = Y_m(i,:)';
obs.update(yk, uk);
seq = [
    0 0 0 0 0 0 0 1
    0 0 0 0 1 0 0 1
    0 0 0 1 0 0 0 1
    0 0 0 0 0 0 1 1
    0 0 0 0 0 1 0 1
];
%disp(obs.i)
%disp(debug_array(obs))
assert(isequaln(obs.i, [7 1]))
assert(isequaln(obs.i_next, [8 1]))
assert(isequaln(cell2mat(obs.seq), seq))
assert(isequaln(obs.f_hold, [4 5]))
assert(isequaln(obs.f_main, [1 2 3]))

% Check probabilities
assert(isequal(obs.gamma_k, [0 0 0 1 0]'))
assert(isequal(obs.p_gamma_k, [0.99 0.99 0.99 0.01 0.99]'))
assert(isequal(round(obs.p_yk_g_seq_Ykm1, 4), [0.1865 0.8015 0.6345 0.1865 0.1865]'))
assert(isequal(round(obs.p_seq_g_Ykm1, 4), [0.9638 0.0097 0.0039 0.0097 0.0097]'))
assert(isequal(round(obs.p_seq_g_Yk, 4), [0.9281 0.0403 0.0129 0.0094 0.0094]'))

% For comparison: probability densities if yk had been 0:
% assert(isequal(round(obs.p_yk_g_seq_Ykm1, 4), [3.1646 1.2260 3.1646 3.1646 3.1646]'))
% assert(isequal(round(obs.p_seq_g_Yk, 4), [0.0050    0.0019    0.0050    0.4940    0.4940]'))

% Check estimates
assert(isequal(round(obs.xkp1_est, 4), [0.4290; 0.1510]))
assert(isequal(round(obs.ykp1_est, 4), 0.1287))
assert(isequal(round(obs.ykp1_est - yk, 4), -0.1713))

% Update at k = 7
i = 8;
uk = U_m(i,:)';
yk = Y_m(i,:)';
obs.update(yk, uk);
seq = [
    0 0 0 0 0 0 0 0
    0 0 0 0 1 0 0 0
    0 0 0 0 0 0 0 1
    0 0 0 0 0 0 1 0
    0 0 0 0 0 1 0 0
];
%disp(obs.i)
%disp(debug_array(obs))
assert(isequaln(obs.i, [8 1]))
assert(isequaln(obs.i_next, [1 1]))
assert(isequaln(cell2mat(obs.seq), seq))
assert(isequaln(obs.f_hold, [3 4]))
assert(isequaln(obs.f_main, [1 2 5]))

% Check probabilities
assert(isequal(obs.gamma_k, [0 0 1 0 0]'))
assert(isequal(obs.p_gamma_k, [0.99 0.99 0.01 0.99 0.99]'))
assert(isequal(round(obs.p_yk_g_seq_Ykm1, 4), [0.0166 1.9388 0.0166 0.0166 0.5807]'))
assert(isequal(round(obs.p_seq_g_Ykm1, 4), [0.9188 0.0399 0.0093 0.0093 0.0093]'))
assert(isequal(round(obs.p_seq_g_Yk, 4), [0.1553 0.7867 0.0016 0.0016 0.0548]'))

% Check estimates
assert(isequal(round(obs.xkp1_est, 4), [1.9184; 0.8599]))
assert(isequal(round(obs.ykp1_est, 4), 0.5755))
assert(isequal(round(obs.ykp1_est - yk, 4), 0.0655))  % old version was 0.0486

% Update at k = 8
i = 9;
uk = U_m(i,:)';
yk = Y_m(i,:)';
obs.update(yk, uk);
seq = [
    0 0 0 0 0 0 0 0
    0 0 0 0 1 0 0 0
    0 0 0 0 0 0 0 1
    0 0 0 0 0 0 1 0
    1 0 0 0 1 0 0 0  % New additions have looped back to 1st position
];
%disp(obs.i)
%disp(debug_array(obs))
assert(isequaln(obs.i, [1 1]))
assert(isequaln(obs.i_next, [2 1]))
assert(isequaln(cell2mat(obs.seq), seq))
assert(isequaln(obs.f_hold, [5 3]))
assert(isequaln(obs.f_main, [1 2 4]))

% Check probabilities
assert(isequal(obs.gamma_k, [0 0 0 0 1]'))
assert(isequal(obs.p_gamma_k, [0.99 0.99 0.99 0.99 0.01]'))
assert(isequal(round(obs.p_yk_g_seq_Ykm1, 4), [0.0084 2.5317 0.0084 0.5437 2.5317]'))
assert(isequal(round(obs.p_seq_g_Ykm1, 4), [0.1537 0.7789 0.0016 0.0016 0.0079]'))
assert(isequal(round(obs.p_seq_g_Yk, 4), [0.0006 0.9889 0.0000 0.0004 0.0100]'))

% Check estimates
assert(isequal(round(obs.xkp1_est, 4), [2.4869; 0.9775]))
assert(isequal(round(obs.ykp1_est, 4), 0.7461))
assert(isequal(round(obs.ykp1_est - yk, 4), 0.0891))  % old version was 0.0797


%% Test sequence updates 2

% Load system and disturbance model from file
sys_rodin_step

% Load observers from file
obs_rodin_step

nT = 6;
x0 = [0; 0];
n_filt = 5;
f = 8;
n_min = 1;  % NOTE: this produces identical results to previous
            % MKFObserverSP and mkf_observer_AFMM with n_min = 2;
            % This is due to the interpretation about whether a
            % hypothesis leaving holding group goes into main group
            % first (as in this version) or can be immediately 
            % eliminated before going to main group (as previously).
obs = MKFObserverSP2(A,B,C,D,Ts,u_meas,P0,epsilon,sigma_wp, ...
    Q0,R,n_filt,f,n_min,label,x0);
assert(isequal(obs.xkp1_est, x0))
assert(isequal(obs.ykp1_est, 0))
assert(obs.d == 1)

% % Generate test simulation data
% nT = 10;
% U_m = zeros(nT+1, sum(u_meas));
% % Add a random shock
% Wp = zeros(nT+1, sum(~u_meas));
% Wp(5, :) = 1;
% % Compute outputs (n0 measurement noise)
% [Y_m, t] = lsim(Gpss, [U_m Wp], t);
% [t U_m Wp Y_m]

% Test simuation data
sim_data = [ ...
         0         0         0         0;
    0.5000         0         0         0;
    1.0000         0         0         0;
    1.5000         0         0         0;
    2.0000         0    1.0000         0;
    2.5000         0         0         0;
    3.0000         0         0    0.3000;
    3.5000         0         0    0.5100;
    4.0000         0         0    0.6570;
    4.5000         0         0    0.7599;
    5.0000         0         0    0.8319];
nT = size(sim_data, 1);
t = sim_data(:, 1);
U_m = sim_data(:, 2);
Wp = sim_data(:, 3);
Y_m = sim_data(:, 4);

% Set marker values on each sequence - for testing only
% these values at the end of the sequences are not used
% by the observer.
for i = 1:n_filt
    obs.seq{i}(8) = i;
end
seq0 = [
    0 0 0 0 0 0 0 1
    0 0 0 0 0 0 0 2
    0 0 0 0 0 0 0 3
    0 0 0 0 0 0 0 4
    0 0 0 0 0 0 0 5
];
%%disp(obs.i)  % use for debugging
%disp(debug_array(obs))
assert(isequaln(obs.i, [0 0]))
assert(isequaln(obs.i_next, [1 1]))
assert(isequaln(cell2mat(obs.seq), seq0))
assert(isequal(obs.n_hold, 1))
assert(isequal(obs.n_main, 4))
assert(isequaln(obs.f_hold, 5))
assert(isequaln(obs.f_main, [1 2 3 4]))

% Check probabilities
assert(isequal(obs.gamma_k, [0 0 0 0 0]'))
assert(isequal(obs.p_gamma_k, [0 0 0 0 0]'))
assert(isequal(obs.p_yk_g_seq_Ykm1, [0 0 0 0 0]'))
assert(isequal(obs.p_seq_g_Ykm1, [0 0 0 0 0]'))
assert(isequal(obs.p_seq_g_Yk, [1 zeros(1, 4)]'))

% Check initialization of filters
assert(isequal(obs.filters{1}.P, obs.P0))
assert(isequal(obs.filters{2}.P, 1e10*eye(2)))
assert(isequal(obs.filters{3}.P, 1e10*eye(2)))
assert(isequal(obs.filters{4}.P, 1e10*eye(2)))
assert(isequal(obs.filters{5}.P, 1e10*eye(2)))

% Check estimates
assert(isequal(obs.xkp1_est, [0; 0]))
assert(isequal(obs.ykp1_est, 0))

% Update at k = 0
i = 1;
uk = U_m(i,:)';
yk = Y_m(i,:)';
obs.update(yk, uk);
seq = [
    0 0 0 0 0 0 0 1
    1 0 0 0 0 0 0 1
    0 0 0 0 0 0 0 3
    0 0 0 0 0 0 0 4
    0 0 0 0 0 0 0 5
];
%disp(obs.i)
%disp(debug_array(obs))
assert(isequaln(obs.i, [1 1]))
assert(isequaln(obs.i_next, [2 1]))
assert(isequaln(cell2mat(obs.seq), seq))
assert(isequaln(obs.f_hold, 2))
assert(isequaln(obs.f_main, [1 5 3 4]))

% Check probabilities
assert(isequal(obs.gamma_k, [0 1 0 0 0]'))
assert(isequal(obs.p_gamma_k, [0.99 0.01 0.99 0.99 0.99]'))
assert(isequal(round(obs.p_yk_g_seq_Ykm1, 4), [0.0420 0.0420 0 0 0]'))
assert(isequal(round(obs.p_seq_g_Ykm1, 4), [0.99 0.01 0 0 0]'))
assert(isequal(round(obs.p_seq_g_Yk, 4), [0.99 0.01 0 0 0]'))

% Check estimates
assert(isequal(obs.xkp1_est, [0; 0]))
assert(isequal(obs.ykp1_est, 0))
assert(isequal(obs.ykp1_est - yk, 0))

% Update at k = 1
i = 2;
uk = U_m(i,:)';
yk = Y_m(i,:)';
obs.update(yk, uk);
seq = [
    0 0 0 0 0 0 0 1
    1 0 0 0 0 0 0 1
    0 0 0 0 0 0 0 3
    0 0 0 0 0 0 0 4
    0 1 0 0 0 0 0 1
];
%disp(obs.i)
%disp(debug_array(obs))
assert(isequaln(obs.i, [2 1]))
assert(isequaln(obs.i_next, [3 1]))
assert(isequaln(cell2mat(obs.seq), seq))
assert(isequaln(obs.f_hold, 5))
assert(isequaln(obs.f_main, [1 2 3 4]))

% Check probabilities
assert(isequal(obs.gamma_k, [0 0 0 0 1]'))
assert(isequal(obs.p_gamma_k, [0.99 0.99 0.99 0.99 0.01]'))
assert(isequal(round(obs.p_yk_g_seq_Ykm1, 4), [0.0420 0.0420 0 0 0.0420]'))
assert(isequal(round(obs.p_seq_g_Ykm1, 4), [0.9801 0.0099 0 0 0.0099]'))  % NOTE: doesn't quite sum to 1
assert(isequal(round(obs.p_seq_g_Yk, 4), [0.9802 0.0099 0 0 0.0099]'))

% Check estimates
assert(isequal(obs.xkp1_est, [0; 0]))
assert(isequal(obs.ykp1_est, 0))
assert(isequal(obs.ykp1_est - yk, 0))

% Update at k = 2
i = 3;
uk = U_m(i,:)';
yk = Y_m(i,:)';
obs.update(yk, uk);
seq = [
    0 0 0 0 0 0 0 1
    1 0 0 0 0 0 0 1
    0 0 1 0 0 0 0 1
    0 0 0 0 0 0 0 4
    0 1 0 0 0 0 0 1
];
%disp(obs.i)
%disp(debug_array(obs))
assert(isequaln(obs.i, [3 1]))
assert(isequaln(obs.i_next, [4 1]))
assert(isequaln(cell2mat(obs.seq), seq))
assert(isequaln(obs.f_hold, 3))
assert(isequaln(obs.f_main, [1 2 5 4]))

% Check probabilities
assert(isequal(obs.gamma_k, [0 0 1 0 0]'))
assert(isequal(obs.p_gamma_k, [0.99 0.99 0.01 0.99 0.99]'))
assert(isequal(round(obs.p_yk_g_seq_Ykm1, 4), [1.8682 1.0834 1.8682 1.8680 1.8682]'))
assert(isequal(round(obs.p_seq_g_Ykm1, 4), [0.9704 0.0098 0.0098 0 0.0098]'))  % NOTE: doesn't quite sum to 1
assert(isequal(round(obs.p_seq_g_Yk, 4), [0.9746 0.0057 0.0098 0 0.0098]'))

% Check estimates
assert(isequal(obs.xkp1_est, [0; 0]))
assert(isequal(obs.ykp1_est, 0))
assert(isequal(obs.ykp1_est - yk, 0))

% Update at k = 3
i = 4;
uk = U_m(i,:)';
yk = Y_m(i,:)';
obs.update(yk, uk);
seq = [
    0 0 0 0 0 0 0 1
    1 0 0 0 0 0 0 1
    0 0 1 0 0 0 0 1
    0 0 0 1 0 0 0 1
    0 1 0 0 0 0 0 1  % all sequences are now splits from #1
];
%disp(obs.i)
%disp(debug_array(obs))
assert(isequaln(obs.i, [4 1]))
assert(isequaln(obs.i_next, [5 1]))
assert(isequaln(cell2mat(obs.seq), seq))
assert(isequaln(obs.f_hold, 4))
assert(isequaln(obs.f_main, [1 2 5 3]))

% Check probabilities
assert(isequal(obs.gamma_k, [0 0 0 1 0]'))
assert(isequal(obs.p_gamma_k, [0.99 0.99 0.99 0.01 0.99]'))
assert(isequal(round(obs.p_yk_g_seq_Ykm1, 4), [2.4676 2.0186 2.4676 2.4676 1.1707]'))
assert(isequal(round(obs.p_seq_g_Ykm1, 4), [0.9649 0.0057 0.0097 0.0097 0.0097]'))
assert(isequal(round(obs.p_seq_g_Yk, 4), [0.9711 0.0047 0.0098 0.0098 0.0047]'))

% Check estimates
assert(isequal(obs.xkp1_est, [0; 0]))
assert(isequal(obs.ykp1_est, 0))
assert(isequal(obs.ykp1_est - yk, 0))

% Update at k = 4
i = 5;
uk = U_m(i,:)';
yk = Y_m(i,:)';
obs.update(yk, uk);
seq = [
    0 0 0 0 0 0 0 1
    0 0 0 0 1 0 0 1  % this seq. has now been replaced
    0 0 1 0 0 0 0 1
    0 0 0 1 0 0 0 1
    0 1 0 0 0 0 0 1 
];
%disp(obs.i)
%disp(debug_array(obs))
assert(isequaln(obs.i, [5 1]))
assert(isequaln(obs.i_next, [6 1]))
assert(isequaln(cell2mat(obs.seq), seq))
assert(isequaln(obs.f_hold, 2))
assert(isequaln(obs.f_main, [1 4 5 3]))

% Check probabilities
assert(isequal(obs.gamma_k, [0 1 0 0 0]'))
assert(isequal(obs.p_gamma_k, [0.99 0.01 0.99 0.99 0.99]'))
assert(isequal(round(obs.p_yk_g_seq_Ykm1, 4), [2.8078 2.8078 1.2019 2.8078 2.0187]'))
assert(isequal(round(obs.p_seq_g_Ykm1, 4), [0.9614 0.0097 0.0097 0.0097 0.0046]'))
assert(isequal(round(obs.p_seq_g_Yk, 4), [0.9728 0.0098 0.0042 0.0098 0.0034]'))

% Check estimates
assert(isequal(obs.xkp1_est, [0; 0]))
assert(isequal(obs.ykp1_est, 0))
assert(isequal(obs.ykp1_est - yk, 0))

% Update at k = 5
i = 6;
uk = U_m(i,:)';
yk = Y_m(i,:)';
obs.update(yk, uk);
seq = [
    0 0 0 0 0 0 0 1
    0 0 0 0 1 0 0 1
    0 0 1 0 0 0 0 1
    0 0 0 1 0 0 0 1
    0 0 0 0 0 1 0 1
];
%disp(obs.i)
%disp(debug_array(obs))
assert(isequaln(obs.i, [6 1]))
assert(isequaln(obs.i_next, [7 1]))
assert(isequaln(cell2mat(obs.seq), seq))
assert(isequaln(obs.f_hold, 5))
assert(isequaln(obs.f_main, [1 4 2 3]))

% Check probabilities
assert(isequal(obs.gamma_k, [0 0 0 0 1]'))
assert(isequal(obs.p_gamma_k, [0.99 0.99 0.99 0.99 0.01]'))
assert(isequal(round(obs.p_yk_g_seq_Ykm1, 4), [3.0218 3.0218 2.0223 1.2172 3.0218]'))
assert(isequal(round(obs.p_seq_g_Ykm1, 4), [0.9631 0.0097 0.0042 0.0097 0.0097]'))
assert(isequal(round(obs.p_seq_g_Yk, 4), [0.9736 0.0098 0.0028 0.0040 0.0098]'))

% Check estimates
assert(isequal(obs.xkp1_est, [0; 0]))
assert(isequal(obs.ykp1_est, 0))
assert(isequal(obs.ykp1_est - yk, 0))

% Update at k = 6  *** First non-zero measurement ***
i = 7;
uk = U_m(i,:)';
yk = Y_m(i,:)';
obs.update(yk, uk);
seq = [
    0 0 0 0 0 0 0 1
    0 0 0 0 1 0 0 1
    0 0 0 0 0 0 1 1
    0 0 0 1 0 0 0 1
    0 0 0 0 0 1 0 1
];
%disp(obs.i)
%disp(debug_array(obs))
assert(isequaln(obs.i, [7 1]))
assert(isequaln(obs.i_next, [8 1]))
assert(isequaln(cell2mat(obs.seq), seq))
assert(isequaln(obs.f_hold, 3))
assert(isequaln(obs.f_main, [1 4 2 5]))

% Check probabilities
assert(isequal(obs.gamma_k, [0 0 1 0 0]'))
assert(isequal(obs.p_gamma_k, [0.99 0.99 0.01 0.99 0.99]'))
assert(isequal(round(obs.p_yk_g_seq_Ykm1, 4), [0.1865 0.8015 0.1865 0.6345 0.1865]'))
assert(isequal(round(obs.p_seq_g_Ykm1, 4), [0.9638 0.0097 0.0097 0.0039 0.0097]'))
assert(isequal(round(obs.p_seq_g_Yk, 4), [0.9281 0.0403 0.0094 0.0129 0.0094]'))

% For comparison: probability densities if yk had been 0:
% assert(isequal(round(obs.p_yk_g_seq_Ykm1, 4), [3.1646 1.2260 3.1646 3.1646 3.1646]'))
% assert(isequal(round(obs.p_seq_g_Yk, 4), [0.0050    0.0019    0.0050    0.4940    0.4940]'))

% Check estimates
assert(isequal(round(obs.xkp1_est, 4), [0.4290; 0.1510]))
assert(isequal(round(obs.ykp1_est, 4), 0.1287))
assert(isequal(round(obs.ykp1_est - yk, 4), -0.1713))

% Update at k = 7
i = 8;
uk = U_m(i,:)';
yk = Y_m(i,:)';
obs.update(yk, uk);
seq = [
    0 0 0 0 0 0 0 0
    0 0 0 0 1 0 0 0
    0 0 0 0 0 0 1 0
    0 0 0 1 0 0 0 0
    0 0 0 0 0 0 0 1
];
%disp(obs.i)
%disp(debug_array(obs))
assert(isequaln(obs.i, [8 1]))
assert(isequaln(obs.i_next, [1 1]))
assert(isequaln(cell2mat(obs.seq), seq))
assert(isequaln(obs.f_hold, 5))
assert(isequaln(obs.f_main, [1 4 2 3]))

% Check probabilities
assert(isequal(obs.gamma_k, [0 0 0 0 1]'))
assert(isequal(obs.p_gamma_k, [0.99 0.99 0.99 0.99 0.01]'))
assert(isequal(round(obs.p_yk_g_seq_Ykm1, 4), [0.0166 1.9388 0.0166 0.9494 0.0166]'))
assert(isequal(round(obs.p_seq_g_Ykm1, 4), [0.9188 0.0399 0.0093 0.0127 0.0093]'))
assert(isequal(round(obs.p_seq_g_Yk, 4), [0.1454 0.7366 0.0015 0.1150 0.0015]'))

% Check estimates
assert(isequal(round(obs.xkp1_est, 4), [1.8619; 0.8147]))
assert(isequal(round(obs.ykp1_est, 4), 0.5586))
assert(isequal(round(obs.ykp1_est - yk, 4), 0.0486))  % same as old version (0.0486)

% Update at k = 8
i = 9;
uk = U_m(i,:)';
yk = Y_m(i,:)';
obs.update(yk, uk);
seq = [
    0 0 0 0 0 0 0 0
    0 0 0 0 1 0 0 0
    1 0 0 0 1 0 0 0
    0 0 0 1 0 0 0 0
    0 0 0 0 0 0 0 1  % New additions have looped back to 1st position
];
%disp(obs.i)
%disp(debug_array(obs))
assert(isequaln(obs.i, [1 1]))
assert(isequaln(obs.i_next, [2 1]))
assert(isequaln(cell2mat(obs.seq), seq))
assert(isequaln(obs.f_hold, 3))
assert(isequaln(obs.f_main, [1 4 2 5]))

% Check probabilities
assert(isequal(obs.gamma_k, [0 0 1 0 0]'))
assert(isequal(obs.p_gamma_k, [0.99 0.99 0.01 0.99 0.99]'))
assert(isequal(round(obs.p_yk_g_seq_Ykm1, 4), [0.0084 2.5317 2.5317 1.3806 0.0084]'))
assert(isequal(round(obs.p_seq_g_Ykm1, 4), [0.1439 0.7293 0.0074 0.1139 0.0015]'))
assert(isequal(round(obs.p_seq_g_Yk, 4), [0.0006 0.9125 0.0092 0.0777 0.0000]'))

% Check estimates
assert(isequal(round(obs.xkp1_est, 4), [2.4556; 0.9601]))
assert(isequal(round(obs.ykp1_est, 4), 0.7367))
assert(isequal(round(obs.ykp1_est - yk, 4), 0.0797))  % Same as old version (0.0797)


%% Run full simulation

clear all
plot_dir = 'plots';

seed = 0;
rng(seed)

% Load system and disturbance model from file
sys_rodin_step

% Load observers from file
obs_rodin_step

% Simulation settings
nT = 100;
t = Ts*(0:nT)';

% Choose time and amplitude of input disturbance
t_shock = 10;
du0 = 1;

% Inputs
%U = (idinput(size(t)) + 1)/2;
U = zeros(size(t));
U(t >= 1) = -1;
% Disturbance input
% This is used by the SKF observer
alpha = zeros(nT+1, 1);
alpha(t == t_shock(1), 1) = 1;
Wp = du0' .* alpha;

% Calculate the input disturbance
% P = zeros(size(U));
% P(t >= t_shock, 1) = du0;

% Combined inputs for simulation
U_sim = [U Wp];

% Simulate system
X = zeros(nT+1,n);
Y = zeros(nT+1,ny);
xk = zeros(n,1);

for i = 1:nT+1

    % Inputs
    uk = U_sim(i, :)';

    % Compute y(k)
    yk = C*xk + D*uk;

    % Store results
    X(i, :) = xk';
    Y(i, :) = yk';
    
    % Compute x(k+1)
    xk = A*xk + B*uk;

end

% Check simulation output is correct
[Y2, t, X2] = lsim(Gpss,U_sim,t);
assert(isequal(X, X2))
assert(isequal(Y, Y2))

% Choose measurement noise for plant
sigma_MP = 0;  % Set to zero for testing
Y_m = Y + sigma_MP'.*randn(size(Y));

% Define custom MKF test observers

% Devise a custom multi-model filter with a shock indicator 
% sequence that perfectly reflects the shock occurence in
% this test simulation (t = t_shock)
% Multiple model filter 1
A2 = repmat({A}, 1, 2);
Bu2 = repmat({Bu}, 1, 2);
C2 = repmat({C}, 1, 2);
Du2 = repmat({Du}, 1, 2);
P0 = 1000*eye(n);
Q0 = diag([q1 0]);
P0_init = repmat({P0}, 1, 2);
Q2 = {diag([Q0(1,1) sigma_wp(1,1)^2]), ...
      diag([Q0(1,1) sigma_wp(1,2)^2])};
R2 = {sigma_M^2, sigma_M^2};
seq = {zeros(1, nT+1); zeros(1, nT+1)};
seq{2}(t == 10) = 1;
p_gamma = [1-epsilon epsilon]';
T = repmat(p_gamma', 2, 1);
d = 1;
MKF3 = MKFObserver(A2,Bu2,C2,Du2,Ts,P0,Q2,R2,seq,T,d,'MKF3');

seq = {zeros(1, nT+1)};
seq{1}(t == 10) = 1;
p_gamma = [1-epsilon epsilon]';
T = repmat(p_gamma', 2, 1);
d = 1;
MKF4 = MKFObserver(A2,Bu2,C2,Du2,Ts,P0,Q2,R2,seq,T,d,'MKF4');

% Define scheduled Kalman filter
% Note: in the case of more than one random input variable, all
% possible combinations of the switching systems need to be 
% accounted for.
% Here, we account for 3 possible combinations:
% combs = [0 0; 1 0; 0 1];
% (This is the same as the MKF filters for the RODD).
% seq = sum(alpha .* 2.^(1:-1:0), 2)';
SKF = MKFObserverSched(A2,Bu2,C2,Du2,Ts,P0,Q2,R2,seq{1},"SKF");

% Simulate observers

% Choose observers to test
observers = {KF2, KF3, SKF, MKF_SP21, MKF3, MKF4};  % , AFMM2, MKF3, MKF4
% Note: KF1 is too slow to pass static error test here

% Measured inputs (not including disturbances)
U_m = U;
n_obs = numel(observers);
MSE = containers.Map();

for i = 1:n_obs

    obs = observers{i};
    [obs, sim_results] = run_test_simulation(nT,Ts,n,ny,U_m,Y_m,obs,alpha);

    % Check observer errors are zero prior to
    % input disturbance
    assert(all(abs(sim_results.X_est(1:20,:) - X(1:20, :)) < 1e-10, [1 2]))
    assert(all(abs(sim_results.Y_est(1:20,:) - Y(1:20, :)) < 1e-10))

    % Check observer static errors are small
    % after input disturbance
    if all(sigma_MP == 0)
        assert(abs(sim_results.Y_est(end, :) - Y(end, :)) < 3e-4);
        assert(abs(sim_results.X_est(end, 2) - du0) < 4e-4);
    end
    % TODO: Errors for AFMM1 were not as low as for RODD MKF observers

    % Compute mean-squared error
    Y_est = sim_results.Y_est;
    MSE(obs.label) = mean((Y_est - Y).^2);
    
    % Save updated observer
    observers{i} = obs;

end

% % Show table of mean-squared errors
% table(MSE.keys', cell2mat(MSE.values'), ...
%     'VariableNames', {'Observer', 'MSE'})

MSE_test_values = containers.Map(...
    {'MMKF_SP21', 'MMKF_SP22', 'KF2', 'KF3', 'SKF', 'MKF3', 'MKF4'}, ...
    [0.002680 0.002687 0.000934 0.003524 0.000929 0.002711 0.000929]' ...
);

for label = MSE.keys
    %fprintf("%s: %f (%f)\n", label{1}, MSE(label{1}), MSE_test_values(label{1}))
    assert(isequal(round(MSE(label{1}), 6), MSE_test_values(label{1})))
end

% Results from old mkf_observer_AFMM code:
% AFMM1: 0.002679 (0.002679)
% KF2: 0.000934 (0.000934)
% KF3: 0.003524 (0.003524)
% SKF: 0.000929 (0.000929)

% % Display results of last simulation
% 
% X_est = sim_results.X_est;
% E_obs = sim_results.E_obs;
% K_obs = sim_results.K_obs;
% trP_obs = sim_results.trP_obs;
% 
% table(t,alpha,U,P,Wp,X,Y,Y_m,X_est,Y_est,E_obs)
% 
% % Display gains and trace of covariance matrix
% table(t, cell2mat(K_obs), cell2mat(trP_obs), ...
%     'VariableNames',{'t', 'K{1}, K{2}', 'trace(P)'})
% 
% % Display AFMM filter groupings
% switch obs.label
%     case {'AFMM1', 'AFMM2'}
%     f_hold = sim_results.AFMM_f_hold
%     f_main = sim_results.AFMM_f_main
%     [array2table(f_hold) array2table(f_main)]
% end
% 
% % Plot of inputs and outputs
% 
% set(groot,'defaultAxesTickLabelInterpreter','latex');
% set(groot,'defaulttextinterpreter','latex');
% set(groot,'defaultLegendInterpreter','latex');
% 
% figure(1); clf
% colors = get(gca,'colororder');
% ax1 = subplot(4,1,1);
% stairs(t,Y_m); hold on
% stairs(t,Y_est,'Linewidth',2);
% ax1.ColorOrder = colors(1:size(Y_m,2),:);
% max_min = [min(min([Y_m Y_est])) max(max([Y_m Y_est]))];
% bd = max([0.1 diff(max_min)*0.1]);
% ylim(max_min + [-bd bd])
% xlabel('t')
% ylabel('$y_m(k)$ and $\hat{y}(k)$')
% title('Process output measurements and estimates')
% legend('$y_m(k)$','$\hat{y}(k)$')
% grid on
% 
% ax2 = subplot(4,1,2);
% stairs(t,X); hold on
% stairs(t,X_est,'Linewidth',2);
% ax2.ColorOrder = colors(size(Y,2)+1:size(Y,2)+size(X,2),:);
% max_min = [min(min([X X_est])) max(max([X X_est]))];
% bd = max([0.1 diff(max_min)*0.1]);
% ylim(max_min + [-bd bd])
% xlabel('t')
% ylabel('$x_i(k)$ and $\hat{x}_i(k)$')
% labels = repmat({''}, 1, n*2);
% for i = 1:n
%     labels{i} = sprintf("$x_{%d}(k)$", i);
% end
% for i = 1:n
%     labels{i+n} = sprintf("$%s{x}_{%d}(k)$", '\hat', i);
% end
% legend(labels)
% title('Actual states and estimates')
% grid on
% 
% ax3 = subplot(4,1,3);
% stairs(t,U,'Linewidth',2); hold on;
% stairs(t,Wp,'Linewidth',2)
% stairs(t,P,'Linewidth',2)
% max_min = [min(min([U Wp P])) max(max([U Wp P]))];
% bd = max([0.1 diff(max_min)*0.1]);
% ylim(max_min + [-bd bd])
% xlabel('t')
% ylabel('$u(k)$, $w_p(k)$ and $p(k)$')
% legend('$u(k)$', '$w_p(k)$', '$p(k)$')
% title('Actual process inputs')
% grid on
% 
% ax4 = subplot(4,1,4);
% stairs(t,alpha,'Color',colors(end,:),'Linewidth',2)
% max_min = [min(min(alpha)) max(max(alpha))];
% bd = max([0.1 diff(max_min)*0.1]);
% ylim(max_min + [-bd bd])
% xlabel('t')
% ylabel('$\gamma(k)$')
% title('Random shock sequence')
% grid on
% 
% linkaxes([ax1, ax2 ax3 ax4], 'x')
% 
% set(gcf,'Position',[100 200 560 600]);
% 
% 
% % Plot of conditional filter probabilities
% switch obs.label
%     case {'MKF1', 'MKF2', 'AFMM1', 'AFMM2'}
%         p_seq_g_Yk = sim_results.MKF_p_seq_g_Yk;
%         % Note: first data points are nans,
%         % ignore last data point to make plot wider
% 
%         figure(11); clf
%         t = Ts*(0:nT)';
%         ax_labels = {'$t$', 'MKF filter ($\Gamma(k)$)', '$Pr(\Gamma(k) \mid Y(k))$'};
%         make_waterfall_plot(t(2:end-1), p_seq_g_Yk(2:end-1, :), [0 1], ...
%             ax_labels, [0 82]);
%         filename = sprintf('rod_mkf_observer_test_pyk_wfplot.pdf');
%         save_fig_to_pdf(fullfile(plot_dir, filename));
%         title('Conditional probabilities of y(k)')
% end
% 
% 
% % Plot of trace of filter covariance matrices
% switch obs.label
%     case {'MKF1', 'MKF2', 'AFMM1', 'AFMM2'}
%         trP_obs = cell2mat(sim_results.trP_obs);
% 
%         figure(12); clf
%         t = Ts*(0:nT)';
%         ax_labels = {'$t$', 'MKF filter ($\Gamma(k)$)', '$Tr(P(k))$'};
%         make_waterfall_plot(t, trP_obs, [0 5], ax_labels, [0 82]);
%         filename = sprintf('rod_mkf_observer_test_trP_wfplot.pdf');
%         save_fig_to_pdf(fullfile(plot_dir, filename));
%         title('Trace of covariance matrices')
% 
% end
% 
% % Plot of filter correction gains (k1)
% switch obs.label
%     case {'MKF1', 'MKF2', 'AFMM1', 'AFMM2'}
%         K_obs = cell2mat(sim_results.K_obs);
%         % Select first gain value onlu
%         K1_obs = K_obs(:,1:2:end);
% 
%         figure(13); clf
%         t = Ts*(0:nT)';
%         ax_labels = {'$t$', 'MKF filter ($\Gamma(k)$)', '$K_1$'};
%         make_waterfall_plot(t, K1_obs, [0 6], ax_labels, [0 82]);
%         filename = sprintf('rod_mkf_observer_test_K_wfplot.png');
%         save_fig_to_pdf(fullfile(plot_dir, filename));
%         title('Filter correction gains (k1)')
% end
% 
% % Plot of final sequence values
% switch obs.label
%     case {'AFMM1', 'AFMM2'}
%         Z = double(cell2mat(obs.seq))';
%         if size(Z, 1) > nT
%             Z = Z(1:nT,:);
%         else
%             Z = [Z(1:obs.i(1),:); Z(obs.i(1)+1:end,:)];
%         end
%         seq_len = size(Z, 1);
%         t = Ts*(nT-seq_len+1:nT)';
% 
%         figure(14); clf
%         ax_labels = {'$t$', 'MKF filter ($\Gamma(k)$)', '$\gamma(k)$'};
%         title('Final filter sequence values')
%         make_waterfall_plot(t,Z,[0 1], ax_labels, [0 82]);
%         filename = sprintf('rod_afmm_filter_test.png');
%         save_fig_to_pdf(fullfile(plot_dir, filename));
% end


%% Test initialization on 2x2 system

% Sample time
Ts = 1;

% Discrete time state space model
A = [ 0.8890       0     1 -0.2;
           0  0.8890  -0.2    1;
           0       0     1    0;
           0       0     0    1];
B = [    1 -0.2  0  0;  % TODO: increase the coupling, -0.5?
      -0.2    1  0  0;
         0    0  1  0;
         0    0  0  1];
C = [ 0.1110 0         0  0;
             0  0.1110 0  0];
D = zeros(2, 4);
Gpss = ss(A,B,C,D,Ts);

% Dimensions
n = size(A, 1);
nu = size(B, 2);
ny = size(C, 1);

% Designate measured input and output signals
u_meas = [true; true; false; false];
y_meas = [true; true];

% Observer model without disturbance noise input
Bu = B(:, u_meas);
Du = D(:, u_meas);
nu = sum(u_meas);
nw = sum(~u_meas);

% Disturbance input (used by SKF observer)
Bw = B(:, ~u_meas);
nw = sum(~u_meas);

% RODD random variable parameters
epsilon = [0.01; 0.01];
sigma_M = [0.1; 0.1];
sigma_wp = [0.01 1; 0.01 1];

% Multiple model observer with sequence pruning 1
label = "MKF_SP1";
P0 = 1000*eye(n);
Q0 = diag([0.01 0.01 0 0]);
R = diag(sigma_M.^2);
f = 10;  % sequence history length
n_filt = 15;  % number of filters
n_min = 5;  % minimum life of cloned filters
MKF_SP21 = MKFObserverSP2(A,B,C,D,Ts,u_meas,P0,epsilon,sigma_wp, ...
    Q0,R,n_filt,f,n_min,label);

% Multiple model observer with sequence pruning 2
label = "MKF_SP2";
P0 = 1000*eye(n);
Q0 = diag([0.01 0.01 0 0]);
R = diag(sigma_M.^2);
f = 10;  % sequence history length
n_filt = 30;  % number of filters
n_min = 10;  % minimum life of cloned filters
MKF_SP22 = MKFObserverSP2(A,B,C,D,Ts,u_meas,P0,epsilon,sigma_wp, ...
    Q0,R,n_filt,f,n_min,label);

% Check observer initialization
assert(isequal(MKF_SP21.epsilon, epsilon))
assert(isequal(MKF_SP21.sigma_wp, sigma_wp))
assert(MKF_SP21.n_filt == 15)
assert(MKF_SP21.n_min == 5)
assert(isequal(MKF_SP21.n_hold, 5*2))
assert(isequal(MKF_SP21.n_main, 5))
assert(isequaln(MKF_SP21.f_hold, 6:15))
assert(isequaln(MKF_SP21.f_main, 1:5))
assert(isequaln(MKF_SP21.i, [0 0]))
assert(MKF_SP21.n == 4)
assert(MKF_SP21.nu == 2)
assert(MKF_SP21.ny == 2)
assert(MKF_SP21.nj == 3)
assert(isequal(MKF_SP21.A{1}, A) && isequal(MKF_SP21.A{2}, A))
assert(isequal(MKF_SP21.B{1}, Bu) && isequal(MKF_SP21.B{2}, Bu))
assert(isequal(MKF_SP21.C{1}, C) && isequal(MKF_SP21.C{2}, C))
assert(isequal(MKF_SP21.D{1}, Du) && isequal(MKF_SP21.D{2}, Du))
assert(MKF_SP21.Ts == Ts)
assert(isequaln(MKF_SP21.u_meas, u_meas))
assert(isequal(MKF_SP21.Q{1}, diag([0.01 0.01 sigma_wp(1, 1)^2 sigma_wp(2, 1)^2])))
assert(isequal(MKF_SP21.Q{2}, diag([0.01 0.01 sigma_wp(1, 2)^2 sigma_wp(2, 1)^2])))
assert(isequal(MKF_SP21.Q{3}, diag([0.01 0.01 sigma_wp(1, 1)^2 sigma_wp(2, 2)^2])))
assert(isequal(MKF_SP21.R{1}, R) && isequal(MKF_SP21.R{2}, R))
assert(numel(MKF_SP21.filters) == MKF_SP21.n_filt)
assert(isequal(size(MKF_SP21.seq), [MKF_SP21.n_filt 1]))
assert(isequal(size(cell2mat(MKF_SP21.seq)), [MKF_SP21.n_filt MKF_SP21.f]))
assert(MKF_SP21.f == size(MKF_SP21.seq{1}, 2))
assert(isequal(size(MKF_SP21.xkp1_est), [n 1]))
assert(isequal(size(MKF_SP21.ykp1_est), [ny 1]))
assert(isequal(round(MKF_SP21.p_gamma, 6), [0.980198; 0.009901; 0.009901]))

% Check observer initialization
assert(isequal(MKF_SP22.epsilon, epsilon))
assert(isequal(MKF_SP22.sigma_wp, sigma_wp))
assert(MKF_SP22.n_filt == 30)
assert(MKF_SP22.n_min == 10)
assert(isequal(MKF_SP22.n_hold, 10*2))
assert(isequal(MKF_SP22.n_main, 10))
assert(isequaln(MKF_SP22.f_hold, 11:30))
assert(isequaln(MKF_SP22.f_main, 1:10))
assert(isequaln(MKF_SP22.i, [0 0]))
assert(MKF_SP22.n == 4)
assert(MKF_SP22.nu == 2)
assert(MKF_SP22.ny == 2)
assert(MKF_SP22.nj == 3)
assert(isequal(MKF_SP22.A{1}, A) && isequal(MKF_SP22.A{2}, A))
assert(isequal(MKF_SP22.B{1}, Bu) && isequal(MKF_SP22.B{2}, Bu))
assert(isequal(MKF_SP22.C{1}, C) && isequal(MKF_SP22.C{2}, C))
assert(isequal(MKF_SP22.D{1}, Du) && isequal(MKF_SP22.D{2}, Du))
assert(MKF_SP22.Ts == Ts)
assert(isequaln(MKF_SP22.u_meas, u_meas))
assert(isequal(MKF_SP22.Q{1}, diag([0.01 0.01 sigma_wp(1, 1)^2 sigma_wp(2, 1)^2])))
assert(isequal(MKF_SP22.Q{2}, diag([0.01 0.01 sigma_wp(1, 2)^2 sigma_wp(2, 1)^2])))
assert(isequal(MKF_SP22.Q{3}, diag([0.01 0.01 sigma_wp(1, 1)^2 sigma_wp(2, 2)^2])))
assert(isequal(MKF_SP22.R{1}, R) && isequal(MKF_SP22.R{2}, R))
assert(numel(MKF_SP22.filters) == MKF_SP22.n_filt)
assert(isequal(size(MKF_SP22.seq), [MKF_SP22.n_filt 1]))
assert(isequal(size(cell2mat(MKF_SP22.seq)), [MKF_SP22.n_filt MKF_SP22.f]))
assert(MKF_SP22.f == size(MKF_SP22.seq{1}, 2))
assert(isequal(size(MKF_SP22.xkp1_est), [n 1]))
assert(isequal(size(MKF_SP22.ykp1_est), [ny 1]))
assert(isequal(round(MKF_SP22.p_gamma, 6), [0.980198; 0.009901; 0.009901]))

% Check optional definition with an initial state estimate works
x0 = [0.1; 0.5; -0.2; -0.4];
MKF_SP_testx0 = MKFObserverSP2(A,B,C,D,Ts,u_meas,P0,epsilon,sigma_wp, ...
    Q0,R,n_filt,f,n_min,label,x0);
assert(isequal(MKF_SP_testx0.xkp1_est, x0))
assert(isequal(MKF_SP_testx0.ykp1_est, C * x0))


%% Test sequence on 2x2 system

% Load system and disturbance model from file
sys_rodin_step_2x2sym

% Load observers from file
obs_rodin_step_2x2

x0 = [0; 0; 0; 0];
f = 8;
n_filt = 10;  % number of filters
n_min = 3;  % minimum life of cloned filters
obs = MKFObserverSP2(A,B,C,D,Ts,u_meas,P0,epsilon,sigma_wp, ...
    Q0,R,n_filt,f,n_min,label,x0);
assert(isequal(obs.xkp1_est, x0))
assert(isequal(obs.ykp1_est, [0; 0]))
assert(obs.d == 1)


% % Sample time
% Ts = 1;
% 
% % Discrete time state space model
% A = [ 0.8890       0     1 -0.2;
%            0  0.8890  -0.2    1;
%            0       0     1    0;
%            0       0     0    1];
% B = [    1 -0.2  0  0;  % TODO: increase the coupling, -0.5?
%       -0.2    1  0  0;
%          0    0  1  0;
%          0    0  0  1];
% C = [ 0.1110 0         0  0;
%              0  0.1110 0  0];
% D = zeros(2, 4);
% Gpss = ss(A,B,C,D,Ts);
% 
% % Dimensions
% n = size(A, 1);
% nu = size(B, 2);
% ny = size(C, 1);
% 
% % Designate measured input and output signals
% u_meas = [true; true; false; false];
% y_meas = [true; true];
% 
% % Observer model without disturbance noise input
% Bu = B(:, u_meas);
% Du = D(:, u_meas);
% nu = sum(u_meas);
% nw = sum(~u_meas);
% 
% % Disturbance input (used by SKF observer)
% Bw = B(:, ~u_meas);
% nw = sum(~u_meas);
% 
% % RODD random variable parameters
% epsilon = [0.01; 0.01];
% sigma_M = [0.1; 0.1];
% sigma_wp = [0.01 1; 0.01 1];


% % Generate test simulation data
% nT = 10;
% t = Ts*(0:nT)';
% U_m = zeros(nT+1, sum(u_meas));
% % Add a random shock
% Wp = zeros(nT+1, sum(~u_meas));
% Wp(5, 1) = 1;
% Wp(8, 1) = 1;
% % Compute outputs (n0 measurement noise)
% [Y_m, t] = lsim(Gpss, [U_m Wp], t);
% [t U_m Wp Y_m]

% Test simuation data
sim_data = [ ...
         0         0         0         0         0         0         0;
    1.0000         0         0         0         0         0         0;
    2.0000         0         0         0         0         0         0;
    3.0000         0         0         0         0         0         0;
    4.0000         0         0    1.0000         0         0         0;
    5.0000         0         0         0         0         0         0;
    6.0000         0         0         0         0    0.1110   -0.0222;
    7.0000         0         0    1.0000         0    0.2097   -0.0419;
    8.0000         0         0         0         0    0.2974   -0.0595;
    9.0000         0         0         0         0    0.4864   -0.0973;
   10.0000         0         0         0         0    0.6544   -0.1309
];
nT = size(sim_data, 1);
t = sim_data(:, 1);
U_m = sim_data(:, 2:3);
Wp = sim_data(:, 4:5);
Y_m = sim_data(:, 6:7);

% Define custom MKF test observers

% Devise a custom multi-model filter with a shock indicator 
% sequence that perfectly reflects the shock occurence in
% this test simulation (t = t_shock)
% Multiple model filter 1
A2 = repmat({A}, 1, 2);
Bu2 = repmat({Bu}, 1, 2);
C2 = repmat({C}, 1, 2);
Du2 = repmat({Du}, 1, 2);
P0 = 1000*eye(n);
Q0 = diag([q1 q1 0 0]);
%P0_init = repmat({P0}, 1, n);
Q2 = {diag([Q0(1,1) Q0(2,2) sigma_wp(:,1)'.^2]), ...
      diag([Q0(1,1) Q0(2,2) sigma_wp(:,2)'.^2])};
R2 = {diag(sigma_M.^2), diag(sigma_M.^2)};
seq = {zeros(1, nT+1); zeros(1, nT+1)};
seq{2}(t == 10) = 1;
p_gamma = [1-epsilon epsilon]';
T = repmat(p_gamma', 2, 1);
d = 1;
MKF3 = MKFObserver(A2,Bu2,C2,Du2,Ts,P0,Q2,R2,seq,T,d,'MKF3');

% Scheduled Kalman filter
% P0 = 1000*eye(n);
% Q0 = diag([Q1 Q2 0 0]);
% R = diag(Radj*sigma_M.^2);
% SKF = KalmanFilter(A,Bu,C,Du,Ts,P0,Q0,R,'SKF');
% SKF.type = 'SKF';
% SKF.Q0 = Q0;
% SKF.Bw = Bw;
% SKF.sigma_wp = sigma_wp;

% Define scheduled Kalman filter
% Note: in the case of more than one random input variable, all
% possible combinations of the switching systems need to be 
% accounted for.
% Here, we account for 3 possible combinations:
% combs = [0 0; 1 0; 0 1];
% (This is the same as the MKF filters for the RODD).
% seq = sum(alpha .* 2.^(1:-1:0), 2)';
SKF = MKFObserverSched(A2,Bu2,C2,Du2,Ts,P0,Q2,R2,seq{1},"SKF");

% Set marker values on each sequence - for testing only
% these values at the end of the sequences are not used
% by the observer.
for i = 1:n_filt
    obs.seq{i}(8) = i;
end
seq0 = [
   0   0   0   0   0   0   0   1
   0   0   0   0   0   0   0   2
   0   0   0   0   0   0   0   3
   0   0   0   0   0   0   0   4
   0   0   0   0   0   0   0   5
   0   0   0   0   0   0   0   6
   0   0   0   0   0   0   0   7
   0   0   0   0   0   0   0   8
   0   0   0   0   0   0   0   9
   0   0   0   0   0   0   0  10
];
%%disp(obs.i)  % use for debugging
%disp(debug_array(obs))
assert(isequaln(obs.i, [0 0]))
assert(isequaln(obs.i_next, [1 1]))
assert(isequaln(cell2mat(obs.seq), seq0))
assert(isequal(obs.n_hold, 6))
assert(isequal(obs.n_main, 4))
assert(isequaln(obs.f_hold, [5 6 7 8 9 10]))
assert(isequaln(obs.f_main, [1 2 3 4]))

% Check probabilities
assert(isequal(obs.gamma_k, [0 0 0 0 0 0 0 0 0 0]'))
assert(isequal(obs.p_gamma_k, [0 0 0 0 0 0 0 0 0 0]'))
assert(isequal(obs.p_yk_g_seq_Ykm1, [0 0 0 0 0 0 0 0 0 0]'))
assert(isequal(obs.p_seq_g_Ykm1, [0 0 0 0 0 0 0 0 0 0]'))
assert(isequal(obs.p_seq_g_Yk, [1 zeros(1, 9)]'))

% Check initialization of filters
assert(isequal(obs.filters{1}.P, obs.P0))
assert(isequal(obs.filters{2}.P, 1e10*eye(4)))
assert(isequal(obs.filters{3}.P, 1e10*eye(4)))
assert(isequal(obs.filters{4}.P, 1e10*eye(4)))
assert(isequal(obs.filters{5}.P, 1e10*eye(4)))

% Check estimates
assert(isequal(obs.xkp1_est, [0; 0; 0; 0]))
assert(isequal(obs.ykp1_est, [0; 0]))

% Update at k = 0
i = 1;
uk = U_m(i,:)';
yk = Y_m(i,:)';
obs.update(yk, uk);
seq = [
   0   0   0   0   0   0   0   1
   1   0   0   0   0   0   0   1
   2   0   0   0   0   0   0   1
   0   0   0   0   0   0   0   4
   0   0   0   0   0   0   0   5
   0   0   0   0   0   0   0   6
   0   0   0   0   0   0   0   7
   0   0   0   0   0   0   0   8
   0   0   0   0   0   0   0   9
   0   0   0   0   0   0   0  10
];
%disp(obs.i)
%disp(debug_array(obs))
assert(isequaln(obs.i, [1 1]))
assert(isequaln(obs.i_next, [2 1]))
assert(isequaln(cell2mat(obs.seq), seq))
assert(isequaln(obs.f_hold, [2 3 5 6 7 8]))
assert(isequaln(obs.f_main, [1 9 10 4]))

% Check probabilities
assert(isequal(obs.gamma_k, [0 1 2 0 0 0 0 0 0 0]'))
assert(isequal(round(obs.p_gamma_k, 4), ...
    [0.9802 0.0099 0.0099 0.9802 0.9802 0.9802 0.9802 0.9802 0.9802 0.9802]'))
assert(isequal(round(obs.p_yk_g_seq_Ykm1, 4), ...
    [0.0129 0.0129 0.0129 0 0 0 0 0 0 0]'))
assert(isequal(round(obs.p_seq_g_Ykm1, 4), ...
    [0.9802 0.0099 0.0099 0 0 0 0 0 0 0]'))
assert(isequal(round(obs.p_seq_g_Yk, 4), ...
    [0.9802 0.0099 0.0099 0 0 0 0 0 0 0]'))

% Check estimates
assert(isequal(obs.xkp1_est, [0; 0; 0; 0]))
assert(isequal(obs.ykp1_est, [0; 0]))
assert(isequal(obs.ykp1_est - yk, [0; 0]))

% Update at k = 1
i = 2;
uk = U_m(i,:)';
yk = Y_m(i,:)';
obs.update(yk, uk);
seq = [
   0   0   0   0   0   0   0   1
   1   0   0   0   0   0   0   1
   2   0   0   0   0   0   0   1
   0   0   0   0   0   0   0   4
   0   0   0   0   0   0   0   5
   0   0   0   0   0   0   0   6
   0   0   0   0   0   0   0   7
   0   0   0   0   0   0   0   8
   0   1   0   0   0   0   0   1
   0   2   0   0   0   0   0   1
];
%disp(obs.i)
%disp(debug_array(obs))
assert(isequaln(obs.i, [2 1]))
assert(isequaln(obs.i_next, [3 1]))
assert(isequaln(cell2mat(obs.seq), seq))
assert(isequaln(obs.f_hold, [9 10 2 3 5 6]))
assert(isequaln(obs.f_main, [1 7 8 4]))

% Check probabilities
assert(isequal(obs.gamma_k, [0 0 0 0 0 0 0 0 1 2]'))
assert(isequal(round(obs.p_gamma_k, 4), ...
    [0.9802 0.9802 0.9802 0.9802 0.9802 0.9802 0.9802 0.9802 0.0099 0.0099]'))
assert(isequal(round(obs.p_yk_g_seq_Ykm1, 4), ...
    [0.0134 0.0134 0.0134 0 0 0 0 0 0.0134 0.0134]'))
assert(isequal(round(obs.p_seq_g_Ykm1, 4), ...
    [0.9608 0.0097 0.0097 0 0 0 0 0 0.0097 0.0097]'))
assert(isequal(round(obs.p_seq_g_Yk, 4), ...
    [0.9612 0.0097 0.0097 0 0 0 0 0 0.0097 0.0097]'))

% Check estimates
assert(isequal(obs.xkp1_est, [0; 0; 0; 0]))
assert(isequal(obs.ykp1_est, [0; 0]))
assert(isequal(obs.ykp1_est - yk, [0; 0]))

% Update at k = 2
i = 3;
uk = U_m(i,:)';
yk = Y_m(i,:)';
obs.update(yk, uk);
seq = [
   0   0   0   0   0   0   0   1
   1   0   0   0   0   0   0   1
   2   0   0   0   0   0   0   1
   0   0   0   0   0   0   0   4
   0   0   0   0   0   0   0   5
   0   0   0   0   0   0   0   6
   0   0   1   0   0   0   0   1
   0   0   2   0   0   0   0   1
   0   1   0   0   0   0   0   1
   0   2   0   0   0   0   0   1
];
%disp(obs.i)
%disp(debug_array(obs))
assert(isequaln(obs.i, [3 1]))
assert(isequaln(obs.i_next, [4 1]))
assert(isequaln(cell2mat(obs.seq), seq))
assert(isequaln(obs.f_hold, [7 8 9 10 2 3]))
assert(isequaln(obs.f_main, [1 5 6 4]))

% Check probabilities
assert(isequal(obs.gamma_k, [0 0 0 0 0 0 1 2 0 0]'))
assert(isequal(round(obs.p_gamma_k, 4), ...
    [0.9802 0.9802 0.9802 0.9802 0.9802 0.9802 0.0099 0.0099 0.9802 0.9802]'))
assert(isequal(round(obs.p_yk_g_seq_Ykm1, 4), ...
    [2.9604 2.6604 2.6604 2.9564 2.9564 2.9564 2.9604 2.9604 2.9604 2.9604]'))
assert(isequal(round(obs.p_seq_g_Ykm1, 4), ...
    [0.9421 0.0095 0.0095 0 0 0 0.0095 0.0095 0.0095 0.0095]'))
assert(isequal(round(obs.p_seq_g_Yk, 4), ...
    [0.9447 0.0086 0.0086 0 0 0 0.0095 0.0095 0.0095 0.0095]'))

% Check estimates
assert(isequal(obs.xkp1_est, [0; 0; 0; 0]))
assert(isequal(obs.ykp1_est, [0; 0]))
assert(isequal(obs.ykp1_est - yk, [0; 0]))

% Update at k = 3
i = 4;
uk = U_m(i,:)';
yk = Y_m(i,:)';
obs.update(yk, uk);
seq = [
   0   0   0   0   0   0   0   1
   1   0   0   0   0   0   0   1
   2   0   0   0   0   0   0   1
   0   0   0   0   0   0   0   4
   0   0   0   1   0   0   0   1
   0   0   0   2   0   0   0   1
   0   0   1   0   0   0   0   1
   0   0   2   0   0   0   0   1
   0   1   0   0   0   0   0   1
   0   2   0   0   0   0   0   1
];
%disp(obs.i)
%disp(debug_array(obs))
assert(isequaln(obs.i, [4 1]))
assert(isequaln(obs.i_next, [5 1]))
assert(isequaln(cell2mat(obs.seq), seq))
assert(isequaln(obs.f_hold, [5 6 7 8 9 10]))
assert(isequaln(obs.f_main, [1 2 3 4]))

% Check probabilities
assert(isequal(obs.gamma_k, [0 0 0 0 1 2 0 0 0 0]'))
assert(isequal(round(obs.p_gamma_k, 4), ...
    [0.9802 0.9802 0.9802 0.9802 0.0099 0.0099 0.9802 0.9802 0.9802 0.9802]'))
assert(isequal(round(obs.p_yk_g_seq_Ykm1, 4), ...
    [5.3037 4.9399 4.9399 5.3019 5.3037 5.3037 5.3037 5.3037 4.4399 4.4399]'))
assert(isequal(round(obs.p_seq_g_Ykm1, 4), ...
    [0.926 0.0084 0.0084 0 0.0094 0.0094 0.0094 0.0094 0.0094 0.0094]'))
assert(isequal(round(obs.p_seq_g_Yk, 4), ...
    [0.9309 0.0079 0.0079 0 0.0094 0.0094 0.0094 0.0094 0.0079 0.0079]'))

% Check estimates
assert(isequal(obs.xkp1_est, [0; 0; 0; 0]))
assert(isequal(obs.ykp1_est, [0; 0]))
assert(isequal(obs.ykp1_est - yk, [0; 0]))

% Update at k = 4
i = 5;
uk = U_m(i,:)';
yk = Y_m(i,:)';
obs.update(yk, uk);
seq = [
   0   0   0   0   0   0   0   1
   0   0   0   0   2   0   0   1  % these seq.s have now been replaced
   2   0   0   0   0   0   0   1  %
   0   0   0   0   1   0   0   1
   0   0   0   1   0   0   0   1
   0   0   0   2   0   0   0   1
   0   0   1   0   0   0   0   1
   0   0   2   0   0   0   0   1
   0   1   0   0   0   0   0   1
   0   2   0   0   0   0   0   1
];
%disp(obs.i)
%disp(debug_array(obs))
assert(isequaln(obs.i, [5 1]))
assert(isequaln(obs.i_next, [6 1]))
assert(isequaln(cell2mat(obs.seq), seq))
assert(isequaln(obs.f_hold, [4 2 5 6 7 8]))
assert(isequaln(obs.f_main, [1 10 3 9]))

% Check probabilities
assert(isequal(obs.gamma_k, [0 2 0 1 0 0 0 0 0 0]'))
assert(isequal(round(obs.p_gamma_k, 4), ...
    [0.9802 0.0099 0.9802 0.0099 0.9802 0.9802 0.9802 0.9802 0.9802 0.9802]'))
assert(isequal(round(obs.p_yk_g_seq_Ykm1, 4), ...
    [7.0443 7.0443 6.7215 7.0443 7.0443 7.0443 5.6272 5.6272 5.9658 5.9658]'))
assert(isequal(round(obs.p_seq_g_Ykm1, 4), ...
    [0.9125 0.0092 0.0077 0.0092 0.0092 0.0092 0.0092 0.0092 0.0077 0.0077]'))
assert(isequal(round(obs.p_seq_g_Yk, 4), ...
    [0.9268 0.0094 0.0075 0.0094 0.0094 0.0094 0.0075 0.0075 0.0066 0.0066]'))

% Check estimates
assert(isequal(obs.xkp1_est, [0; 0; 0; 0]))
assert(isequal(obs.ykp1_est, [0; 0]))
assert(isequal(obs.ykp1_est - yk, [0; 0]))

% Update at k = 5
i = 6;
uk = U_m(i,:)';
yk = Y_m(i,:)';
obs.update(yk, uk);
seq = [
   0   0   0   0   0   0   0   1
   0   0   0   0   2   0   0   1
   2   0   0   0   0   0   0   1
   0   0   0   0   1   0   0   1
   0   0   0   1   0   0   0   1
   0   0   0   2   0   0   0   1
   0   0   1   0   0   0   0   1
   0   0   2   0   0   0   0   1
   0   0   0   0   0   2   0   1
   0   0   0   0   0   1   0   1
];
%disp(obs.i)
%disp(debug_array(obs))
assert(isequaln(obs.i, [6 1]))
assert(isequaln(obs.i_next, [7 1]))
assert(isequaln(cell2mat(obs.seq), seq))
assert(isequaln(obs.f_hold, [10 9 4 2 5 6]))
assert(isequaln(obs.f_main, [1 7 3 8]))

% Check probabilities
assert(isequal(obs.gamma_k, [0 0 0 0 0 0 0 0 2 1]'))
assert(isequal(round(obs.p_gamma_k, 4), ...
    [0.9802 0.9802 0.9802 0.9802 0.9802 0.9802 0.9802 0.9802 0.0099 0.0099]'))
assert(isequal(round(obs.p_yk_g_seq_Ykm1, 4), ...
    [8.3546 8.3546 8.0866 8.3546 6.4601 6.4601 6.5830 6.5830 8.3546 8.3546]'))
assert(isequal(round(obs.p_seq_g_Ykm1, 4), ...
    [0.9085 0.0092 0.0073 0.0092 0.0092 0.0092 0.0073 0.0073 0.0092 0.0092]'))
assert(isequal(round(obs.p_seq_g_Yk, 4), ...
    [0.9289 0.0094 0.0073 0.0094 0.0073 0.0073 0.0059 0.0059 0.0094 0.0094]'))

% Check estimates
assert(isequal(obs.xkp1_est, [0; 0; 0; 0]))
assert(isequal(obs.ykp1_est, [0; 0]))
assert(isequal(obs.ykp1_est - yk, [0; 0]))

% Update at k = 6  *** First non-zero measurement ***
i = 7;
uk = U_m(i,:)';
yk = Y_m(i,:)';
obs.update(yk, uk);
seq = [
   0   0   0   0   0   0   0   1
   0   0   0   0   2   0   0   1
   2   0   0   0   0   0   0   1
   0   0   0   0   1   0   0   1
   0   0   0   1   0   0   0   1
   0   0   0   2   0   0   0   1
   0   0   0   0   0   0   2   1
   0   0   0   0   0   0   1   1
   0   0   0   0   0   2   0   1
   0   0   0   0   0   1   0   1
];
%disp(obs.i)
%disp(debug_array(obs))
assert(isequaln(obs.i, [7 1]))
assert(isequaln(obs.i_next, [8 1]))
assert(isequaln(cell2mat(obs.seq), seq))
assert(isequaln(obs.f_hold, [8 7 10 9 4 2]))
assert(isequaln(obs.f_main, [1 6 3 5]))

% Check probabilities
assert(isequal(obs.gamma_k, [0 0 0 0 0 0 2 1 0 0]'))
assert(isequal(round(obs.p_gamma_k, 4), ...
    [0.9802 0.9802 0.9802 0.9802 0.9802 0.9802 0.0099 0.0099 0.9802 0.9802]'))
assert(isequal(round(obs.p_yk_g_seq_Ykm1, 4), ...
    [6.4234 4.9679 6.2889 5.7034 5.6712 4.9222 6.4234 6.4234 6.4234 6.4234]'))
assert(isequal(round(obs.p_seq_g_Ykm1, 4), ...
    [0.9105 0.0092 0.0071 0.0092 0.0071 0.0071 0.0092 0.0092 0.0092 0.0092]'))
assert(isequal(round(obs.p_seq_g_Yk, 4), ...
    [0.9279 0.0072 0.0071 0.0083 0.0064 0.0056 0.0094 0.0094 0.0094 0.0094]'))

% Check estimates
assert(isequal(round(obs.xkp1_est, 4), [0.4829; -0.0994; 0.1130; -0.0017]))
assert(isequal(round(obs.ykp1_est, 4), [0.0536; -0.0110]))
assert(isequal(round(obs.ykp1_est - yk, 4), [-0.0574; 0.0112]))

% Update at k = 7
i = 8;
uk = U_m(i,:)';
yk = Y_m(i,:)';
obs.update(yk, uk);
seq = [
   0   0   0   0   0   0   0   0
   0   0   0   0   2   0   0   0
   2   0   0   0   0   0   0   0
   0   0   0   0   1   0   0   0
   0   0   0   0   0   0   0   2
   0   0   0   0   0   0   0   1
   0   0   0   0   0   0   2   0
   0   0   0   0   0   0   1   0
   0   0   0   0   0   2   0   0
   0   0   0   0   0   1   0   0
];
%disp(obs.i)
%disp(debug_array(obs))
assert(isequaln(obs.i, [8 1]))
assert(isequaln(obs.i_next, [1 1]))
assert(isequaln(cell2mat(obs.seq), seq))
assert(isequaln(obs.f_hold, [6 5 8 7 10 9]))
assert(isequaln(obs.f_main, [1 4 3 2]))

% Check probabilities
assert(isequal(obs.gamma_k, [0 0 0 0 2 1 0 0 0 0]'))
assert(isequal(round(obs.p_gamma_k, 4), ...
    [0.9802 0.9802 0.9802 0.9802 0.0099 0.0099 0.9802 0.9802 0.9802 0.9802]'))
assert(isequal(round(obs.p_yk_g_seq_Ykm1, 4), ...
    [4.4809 3.5632 4.4484 6.3667 4.4809 4.4809 4.4809 4.4809 3.509 4.803]'))
assert(isequal(round(obs.p_seq_g_Ykm1, 4), ...
    [0.9095 0.0071 0.007 0.0082 0.0092 0.0092 0.0092 0.0092 0.0092 0.0092]'))
assert(isequal(round(obs.p_seq_g_Yk, 4), ...
    [0.9211 0.0057 0.007 0.0117 0.0093 0.0093 0.0093 0.0093 0.0073 0.01]'))

% Check estimates
assert(isequal(round(obs.xkp1_est, 4), [1.1284; -0.2311; 0.2404; -0.0031]))
assert(isequal(round(obs.ykp1_est, 4), [0.1253; -0.0256]))
assert(isequal(round(obs.ykp1_est - yk, 4), [-0.0844; 0.0163]))  % same as old version (0.0486)

% Update at k = 8
i = 9;
uk = U_m(i,:)';
yk = Y_m(i,:)';
obs.update(yk, uk);
seq = [
   0   0   0   0   0   0   0   0
   1   0   0   0   0   0   0   0
   2   0   0   0   0   0   0   0
   0   0   0   0   1   0   0   0
   0   0   0   0   0   0   0   2
   0   0   0   0   0   0   0   1
   0   0   0   0   0   0   2   0
   0   0   0   0   0   0   1   0
   0   0   0   0   0   2   0   0
   0   0   0   0   0   1   0   0  % New additions have looped back to 1st position
];
%disp(obs.i)
%disp(debug_array(obs))
assert(isequaln(obs.i, [1 1]))
assert(isequaln(obs.i_next, [2 1]))
assert(isequaln(cell2mat(obs.seq), seq))
assert(isequaln(obs.f_hold, [2 3 6 5 8 7]))
assert(isequaln(obs.f_main, [1 4 9 10]))

% Check probabilities
assert(isequal(obs.gamma_k, [0 1 2 0 0 0 0 0 0 0]'))
assert(isequal(round(obs.p_gamma_k, 4), ...
    [0.9802 0.0099 0.0099 0.9802 0.9802 0.9802 0.9802 0.9802 0.9802 0.9802]'))
assert(isequal(round(obs.p_yk_g_seq_Ykm1, 4), ...
    [3.677 3.677 3.677 7.9723 3.677 3.677 2.896 4.4379 2.9727 6.9293]'))
assert(isequal(round(obs.p_seq_g_Ykm1, 4), ...
    [0.9028 0.0091 0.0091 0.0115 0.0091 0.0091 0.0091 0.0091 0.0071 0.0098]'))
assert(isequal(round(obs.p_seq_g_Yk, 4), ...
    [0.8969 0.0091 0.0091 0.0248 0.0091 0.0091 0.0071 0.0109 0.0057 0.0183]'))

% Check estimates
assert(isequal(round(obs.xkp1_est, 4), [1.8414; -0.3752; 0.3694; -0.0038]))
assert(isequal(round(obs.ykp1_est, 4), [0.2044; -0.0416]))
assert(isequal(round(obs.ykp1_est - yk, 4), [-0.0930; 0.0179]))  % Same as old version (0.0797)


%% Full simulation on 2x2 system

% Sample time
Ts = 1;

% Discrete time state space model
A = [ 0.8890       0     1 -0.2;
           0  0.8890  -0.2    1;
           0       0     1    0;
           0       0     0    1];
B = [    1 -0.2  0  0;  % TODO: increase the coupling, -0.5?
      -0.2    1  0  0;
         0    0  1  0;
         0    0  0  1];
C = [ 0.1110 0         0  0;
             0  0.1110 0  0];
D = zeros(2, 4);
Gpss = ss(A,B,C,D,Ts);

% Dimensions
n = size(A, 1);
nu = size(B, 2);
ny = size(C, 1);

% Designate measured input and output signals
u_meas = [true; true; false; false];
y_meas = [true; true];

% Observer model without disturbance noise input
Bu = B(:, u_meas);
Du = D(:, u_meas);
nu = sum(u_meas);
nw = sum(~u_meas);

% Disturbance input (used by SKF observer)
Bw = B(:, ~u_meas);
nw = sum(~u_meas);

% RODD random variable parameters
epsilon = [0.01; 0.01];
sigma_M = [0.1; 0.1];
sigma_wp = [0.01 1; 0.01 1];

% Kalman filter 3 - manually tuned
% Covariance matrices
P0 = 1000*eye(n);
Q = diag([0.01 0.01 0.1^2 0.1^2]);
R = diag(sigma_M.^2);
KF3 = KalmanFilter(A,Bu,C,Du,Ts,P0,Q,R,'KF3');

% Multiple model observer with sequence pruning 1
label = "MKF_SP1";
P0 = 1000*eye(n);
Q0 = diag([0.01 0.01 0 0]);
R = diag(sigma_M.^2);
f = 10;  % sequence history length
n_filt = 15;  % number of filters
n_min = 5;  % minimum life of cloned filters
MKF_SP21 = MKFObserverSP2(A,B,C,D,Ts,u_meas,P0,epsilon,sigma_wp, ...
    Q0,R,n_filt,f,n_min,label);

% Multiple model observer with sequence pruning 2
label = "MKF_SP2";
P0 = 1000*eye(n);
Q0 = diag([0.01 0.01 0 0]);
R = diag(sigma_M.^2);
f = 10;  % sequence history length
n_filt = 30;  % number of filters
n_min = 10;  % minimum life of cloned filters
MKF_SP22 = MKFObserverSP2(A,B,C,D,Ts,u_meas,P0,epsilon,sigma_wp, ...
    Q0,R,n_filt,f,n_min,label);

% Simulation settings
nT = 200;
t = Ts*(0:nT)';

% Choose time and amplitude of input disturbance
t_shock = [5 10];
du0 = [1; 1];
% When you make the shock larger the MKF observers
% do better
%du0 = [2; 2];

% Measured input
%U = (idinput(size(t)) + 1)/2;
U = zeros(nT+1, 2);
U(t >= 1, 1) = -1;

% Disturbance input
% This is used by the SKF observer
alpha = zeros(nT+1, 2);
alpha(t == t_shock(1), 1) = 1;
alpha(t == t_shock(2), 2) = 1;
Wp = du0' .* alpha;

U_sim = [U Wp];

% Custom MKF test observer
% Devise a custom multi-model filter with a shock indicator 
% sequence that perfectly reflects the shock occurence in
% this test simulation (t = t_shock)
% Multiple model filter 1
A2 = repmat({A}, 1, 3);
Bu2 = repmat({Bu}, 1, 3);
C2 = repmat({C}, 1, 3);
Du2 = repmat({Du}, 1, 3);
P0 = 1000*eye(n);
Q0 = diag([0.01 0.01 1 1]);
%P0_init = repmat({P0}, 1, 3);
Q2 = {diag([Q0(1,1) Q0(2,2) sigma_wp(1,1)^2 sigma_wp(2,1)^2]), ...
      diag([Q0(1,1) Q0(2,2) sigma_wp(1,2)^2 sigma_wp(2,1)^2]), ...
      diag([Q0(1,1) Q0(2,2) sigma_wp(1,1)^2 sigma_wp(2,2)^2])};
R2 = {diag(sigma_M.^2), diag(sigma_M.^2), diag(sigma_M.^2)};
seq = {zeros(1, nT+1); zeros(1, nT+1); zeros(1, nT+1)};
seq{2}(t == t_shock(1)) = 1;
seq{3}(t == t_shock(2)) = 2;
p_gamma = [1-epsilon epsilon]';
Z = [0 0; 0 1; 1 0];  % combinations
p_gamma = prod(prob_gamma(Z', p_gamma), 1)';
p_gamma = p_gamma ./ sum(p_gamma);  % normalized
T = repmat(p_gamma', 3, 1);
d = 1;
MKF3 = MKFObserver(A2,Bu2,C2,Du2,Ts,P0,Q2,R2,seq,T,d,'MKF3');

seq = {zeros(1, nT+1)};
seq{1}(t == t_shock(1)) = 1;
seq{1}(t == t_shock(2)) = 2;
MKF4 = MKFObserver(A2,Bu2,C2,Du2,Ts,P0,Q2,R2,seq,T,d,'MKF4');

% Define scheduled Kalman filter
% Note: in the case of more than one random input variable, all
% possible combinations of the switching systems need to be 
% accounted for.
% Here, we account for 3 possible combinations:
% combs = [0 0; 1 0; 0 1];
% (This is the same as the MKF filters for the RODD).
% seq = sum(alpha .* 2.^(1:-1:0), 2)';
SKF = MKFObserverSched(A2,Bu2,C2,Du2,Ts,P0,Q2,R2,seq{1},"SKF");

% Choose observers to test
observers = {KF3, MKF_SP21, MKF_SP22, SKF};

% Simulate system
X = zeros(nT+1,n);
Y = zeros(nT+1,ny);
xk = zeros(n,1);

for i = 1:nT+1

    % Inputs
    uk = U_sim(i,:)';

    % Compute y(k)
    yk = C*xk + D*uk;

    % Store results
    X(i, :) = xk';
    Y(i, :) = yk';
    
    % Compute x(k+1)
    xk = A*xk + B*uk;

end

% Check simulation output is correct
[Y2, t, X2] = lsim(Gpss, U_sim, t);
assert(isequal(X, X2))
assert(isequal(Y, Y2))

% Choose measurement noise for plant
sigma_MP = [0; 0];  % Set to zero for testing
Y_m = Y + sigma_MP'.*randn(nT+1, ny);

% Simulate observers

% Measured inputs (not including disturbances)
U_m = U;

n_obs = numel(observers);
MSE = containers.Map();
for i = 1:n_obs

    obs = observers{i};
    [obs, sim_results] = run_test_simulation(nT,Ts,n,ny,U_m,Y_m,obs,alpha);

    % Check observer errors are zero prior to
    % input disturbance
    assert(all(abs(sim_results.X_est(1:5,:) - X(1:5, :)) < 1e-10, [1 2]))
    assert(all(abs(sim_results.Y_est(1:5,:) - Y(1:5, :)) < 1e-10, [1 2]))

    % Check observer static errors are small
    % after input disturbance
    % TODO: Should these be closer?
    if all(sigma_MP == 0)
        assert(all(abs(sim_results.Y_est(end, :) - Y(end, :)) < 1e-3, [1 2]));
        assert(all(abs(sim_results.X_est(end, 3:4) - du0) < 1e-3, [1 2]));
    end

    % Compute mean-squared error
    Y_est = sim_results.Y_est;
    MSE(obs.label) = mean((Y_est - Y).^2);
    %fprintf("%d, %s: %f\n", i, obs.label, mean((Y_est - Y).^2))

    % Save updated observer
    observers{i} = obs;

end


% % Display results of last simulation
% 
% X_est = sim_results.X_est;
% E_obs = sim_results.E_obs;
% K_obs = sim_results.K_obs;
% trP_obs = sim_results.trP_obs;
% 
% table(t,alpha,U,Wp,X,Y,Y_m,X_est,Y_est,E_obs)
% 
% % Display gains and trace of covariance matrix
% table(t, cell2mat(K_obs), cell2mat(trP_obs), ...
%     'VariableNames', {'t', 'K{1}, K{2}', 'trace(P{1}), trace(P{2})'})
% 
% % Show table of mean-squared errors
% table(MSE.keys', cell2mat(MSE.values'), ...
%     'VariableNames', {'Observer', 'MSE'})

% Results on Nov 8 after reverting back the Bayesian updating
MSE_test_values = containers.Map(...
 {'KF3',               'MKF_SP1',              'MKF_SP2',              ...
  'SKF'}, ...
 {[0.000676 0.000936], [0.000739 0.001498], [0.000759 0.001521], ...
  [0.000123 0.000132]} ...
);

for label = MSE.keys
    %fprintf("%s: %f, %f (%f, %f)\n", label{1}, MSE(label{1}), MSE_test_values(label{1}))
    assert(isequal(round(MSE(label{1}), 6), MSE_test_values(label{1})))
end

return


function [obs, sim_results] = run_test_simulation(nT,Ts,n,ny,U_m,Y_m, ...
    obs,alpha)

    k = (0:nT)';
    t = Ts*k;
    X_est = nan(nT+1,n);
    Y_est = nan(nT+1,ny);
    E_obs = nan(nT+1,ny);

    % Arrays to store observer variables
    switch obs.type
        case {'MKF', 'MKF-SF'}
            n_filt = obs.n_filt;
            MKF_p_seq_g_Yk = nan(nT+1, n_filt);
        case {'MKF-SP', 'MKF-SP'}
            n_filt = obs.n_filt;
            MKF_p_seq_g_Yk = nan(nT+1, n_filt);
            AFMM_f_main = nan(nT+1, numel(obs.f_main));
            AFMM_f_hold = nan(nT+1, numel(obs.f_hold));
        otherwise
            n_filt = 1;
    end
    K_obs = cell(nT+1, n_filt);
    trP_obs = cell(nT+1, n_filt);

    % Start simulation at k = 0
    for i = 1:nT+1

        % For debugging:
        %fprintf("t = %f\n", t(i));

        % Process measurements
        uk_m = U_m(i,:)';
        yk_m = Y_m(i,:)';

        % Record observer estimates and output errors
        X_est(i, :) = obs.xkp1_est';
        Y_est(i, :) = obs.ykp1_est';
        E_obs(i, :) = yk_m' - obs.ykp1_est';

        % Kalman update equations
        % Update observer gains and covariance matrix
        switch obs.type

            case {'KF', 'SKF'}
                obs.update(yk_m, uk_m);

                % Record filter gain and covariance matrix
                K_obs{i, 1} = obs.K';
                trP_obs{i, 1} = trace(obs.P);

            case {'MKF', 'MKF-SF'}
                obs.update(yk_m, uk_m);

                % Record filter gains and covariance matrices
                for j=1:obs.n_filt
                    K_obs{i, j} = obs.filters{j}.K';
                    trP_obs{i, j} = trace(obs.filters{j}.P);
                end

                % Record filter conditional probabilities
                MKF_p_seq_g_Yk(i, :) = obs.p_seq_g_Yk';

            case {'MKF-SP'}
                obs.update(yk_m, uk_m);

                % Record filter gains and covariance matrices
                for j=1:obs.n_filt
                    K_obs{i, j} = obs.filters{j}.K';
                    trP_obs{i, j} = trace(obs.filters{j}.P);
                end

                % Record filter conditional probabilities
                MKF_p_seq_g_Yk(i, :) = obs.p_seq_g_Yk';

                % Record filter arrangement
                AFMM_f_main(i, :) = obs.f_main;
                AFMM_f_hold(i, :) = obs.f_hold;

            otherwise
                error('Observer type not valid')

        end

    end

    sim_results.t = t;
    sim_results.k = k;
    sim_results.X_est = X_est;
    sim_results.Y_est = Y_est;
    sim_results.E_obs = E_obs;
    sim_results.K_obs = K_obs;
    sim_results.trP_obs = trP_obs;
    switch obs.type
        case {'MKF', 'MKF-SF', 'MKF-SP'}
            sim_results.MKF_p_seq_g_Yk = MKF_p_seq_g_Yk;
    end
    switch obs.type
        case 'MKF-SP'
            sim_results.AFMM_f_main = AFMM_f_main;
            sim_results.AFMM_f_hold = AFMM_f_hold;
    end

end


function dba = debug_array(obs)
% For debugging and testing sequences
    hold = zeros(obs.n_filt, 1);
    hold(nonzeros(obs.f_hold)) = nonzeros(obs.f_hold);
    main = zeros(obs.n_filt, 1);
    main(nonzeros(obs.f_main)) = nonzeros(obs.f_main);
    seq = cell2mat(obs.seq);
    p_max = (obs.p_seq_g_Yk == max(obs.p_seq_g_Yk));
    dba = [table(hold, main) array2table(seq) table(p_max)];
end