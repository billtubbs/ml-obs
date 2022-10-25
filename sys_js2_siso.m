% Define simple switching system - used for testing

% Sample period
Ts = 0.5;

% Discrete time state space models
% Model #1
model1.A = 0.7;
model1.B = 1;
model1.C = 0.3;
model1.Ts = Ts;
Gpss1 = ss(model1.A,model1.B,model1.C,0,Ts);

% Model #2
model2.A = 0.9;
model2.B = 1;
model2.C = -0.3;  % -ve gain!
model2.Ts = Ts;
Gpss2 = ss(model2.A,model2.B,model2.C,0,Ts);

% Array of system models
models = {model1, model2};

% Dimensions
[nj, n, nu, ny, Ts] = check_models(models);
