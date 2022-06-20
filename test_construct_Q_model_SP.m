% Test construct_Q_model_SP.m function
% (used by MKFObserverSP)

clear all


%% Test on SISO system
% with one randomly-occuring step input disturbance

Q0 = [ 0.0100         0
            0         0];
Bw = [0; 1];
epsilon = 0.01;
sigma_wp = [    0.0100    1.0000];
var_wp = sigma_wp.^2;
nw = 1;

[Q, p_gamma] = construct_Q_model_SP(Q0, Bw, epsilon, var_wp, nw);
Q_test = {
    [ 0.0100       0
           0       0.0001] ... 
    [ 0.0100       0
           0       1     ]
};

assert(isequal(Q, Q_test))
assert(isequal(p_gamma, [0.99; 0.01]))


%% Test on 2x2 system
% with two randomly-occuring step input disturbances

Q0 = [
    0.0100         0         0         0
         0    0.0100         0         0
         0         0         0         0
         0         0         0         0
];

Bw = [
     0     0
     0     0
     1     0
     0     1
];
epsilon = [0.01; 0.01];
sigma_wp = [
    0.0100    1.0000
    0.0100    1.0000
];

d = 1;
var_wp = sigma_wp.^2 / d;
alpha = (1 - (1 - epsilon).^d);
nw = size(sigma_wp, 1);

[Q, p_gamma] = construct_Q_model_SP(Q0, Bw, alpha, var_wp, nw);
Q_test = {
    [    0.0100         0         0         0
              0    0.0100         0         0
              0         0    0.0001         0
              0         0         0    0.0001] ...
    [    0.0100         0         0         0
              0    0.0100         0         0
              0         0    1.0000         0
              0         0         0    0.0001] ...
    [    0.0100         0         0         0
              0    0.0100         0         0
              0         0    0.0001         0
              0         0         0    1.0000]
};
assert(isequal(Q, Q_test));


%% Test on 3 input system
% with three randomly-occuring step input disturbances

Q0 = [
    0.0100         0         0         0         0
         0    0.0100         0         0         0
         0         0         0         0         0
         0         0         0         0         0
         0         0         0         0         0
];

Bw = [
     0     0     0
     0     0     0
     1     0     0
     0     1     0
     0     0     1
];
epsilon = [0.005; 0.0025; 0.0025];
sigma_wp = [
    0.0100    1.0000
    0.0050    0.5000
    0.0050    0.5000
];

d = 1;
var_wp = sigma_wp.^2 / d;
alpha = (1 - (1 - epsilon).^d);
nw = size(sigma_wp, 1);

[Q, p_gamma] = construct_Q_model_SP(Q0, Bw, alpha, var_wp, nw);
Q_test = {
    [    0.0100         0         0         0         0
              0    0.0100         0         0         0
              0         0    0.0001         0         0
              0         0         0    2.5e-5         0
              0         0         0         0    2.5e-5] ...
    [    0.0100         0         0         0         0
              0    0.0100         0         0         0
              0         0         1         0         0
              0         0         0    2.5e-5         0
              0         0         0         0    2.5e-5] ...
    [    0.0100         0         0         0         0
              0    0.0100         0         0         0
              0         0    0.0001         0         0
              0         0         0      0.25         0
              0         0         0         0    2.5e-5] ...
    [    0.0100         0         0         0         0
              0    0.0100         0         0         0
              0         0    0.0001         0         0
              0         0         0    2.5e-5         0
              0         0         0         0      0.25]
};
assert(isequal(Q, Q_test));

p_gamma_test = [
    prod(1-epsilon);
    epsilon(1)*prod(1-epsilon(2:3));
    epsilon(2)*prod(1-epsilon([1, 3]));
    epsilon(3)*prod(1-epsilon(1:2));
];
p_gamma_test = p_gamma_test / sum(p_gamma_test);
assert(max(abs(p_gamma - p_gamma_test)) < 1e-15)
