% Test get_MPC_KF_params_OD.m and get_MPC_KF_params_ID.m


%% Test get_MPC_KF_params_OD on SISO example

% Example from documentation
% https://www.mathworks.com/help/mpc/ug/design-estimator-equivalent-to-mpc-built-in-kf.html
x0 = [1.6931; 4.3863];
u0 = 7.7738;
y0 = 0.5;
G = ss([-2 1;-9.7726 -1.3863],[0;1],[0.5 0],0);
Ts = 0.1;
Gd = c2d(G,Ts);

% Turn off information display
old_status = mpcverbosity('off');

% Define default MPC
mpcobj = mpc(Gd,Ts);

% Set some parameters
mpcobj.Model.Nominal.Y = y0;
mpcobj.Model.Nominal.U = u0;
mpcobj.Model.Nominal.X = x0;
mpcobj.Weights.OutputVariables = 0.2;

% Calculate estimator parameters
[Q, R, Gpred] = get_MPC_KF_params_OD(mpcobj);
assert(isequal(round(Q, 6), [ ...
    0.000020    0.000408       0
    0.000408    0.008453       0
         0         0        0.01
]))
assert(R == 1)
assert(isequal(round(Gpred.A, 6), [ ...
    0.778227    0.083069         0
   -0.811801    0.829206         0
           0           0         1
]))
assert(isequal(round(Gpred.B, 6), [ ...
    0.004435
    0.091939
           0
]))
assert(isequal(round(Gpred.C, 6), [ ...
    0.5  0   1
]))

% Repeat with MPC which has a continuous-time model
mpcobj = mpc(G,Ts);
mpcobj.Model.Nominal.Y = y0;
mpcobj.Model.Nominal.U = u0;
mpcobj.Model.Nominal.X = x0;
mpcobj.Weights.OutputVariables = 0.2;

% Calculate estimator parameters
[Q, R, Gpred] = get_MPC_KF_params_OD(mpcobj);
assert(isequal(round(Q, 6), [ ...
    0.000020    0.000408         0
    0.000408    0.008453         0
         0         0        0.0100
]))
assert(R == 1)
assert(isequal(round(Gpred.A, 6), [ ...
    0.778227    0.083069         0
   -0.811801    0.829206         0
           0           0         1
]))
assert(isequal(round(Gpred.B, 6), [ ...
    0.004435
    0.091939
           0
]))
assert(isequal(round(Gpred.C, 6), [ ...
    0.5  0   1
]))

% Return information display settings to previous value
mpcverbosity(old_status);


%% Test get_MPC_KF_params_OD on 2x2 system

s=tf('s');
G11=2.9*exp(-2*s)/(1+24*s);
G12=2.9*exp(-4*s)/(1+12*s);
G21=-1.74*exp(-6*s)/(1+14*s);
G22=1.74*exp(-4*s)/(1+5*s);
G = [G11 G12; G21 G22];
Gcont=ss(G);
Ts = 2; % Sampling period        
model = c2d(Gcont,Ts);
model = absorbDelay(model);

mpcverbosity off;
model = setmpcsignals(model,'MV',1:2); % inputs 1 and 2 are MVs
model = setmpcsignals(model,'MO',1:2); % outputs 1 and 2 are MOs 
mpcobj = mpc(model,Ts);

% Operating points
u_op = [10;20];
y_op = [30;40];
mpcobj.model.Nominal.U=u_op;
mpcobj.model.Nominal.Y=y_op;

% Prediction and control horizons
mpcobj.PredictionHorizon = 20;
mpcobj.ControlHorizon = 1;

% TODO: Need to get it working with scaling factors on MVs. 'ScaleFactor',{50,50}

% Constraints and scaling for Mvs
mpcobj.MV = struct('Min',{0,0},'Max',{50,50}, ... % u_1min u_2min   u1max u2max
                   'RateMin',{-10,-10},'RateMax',{10,10}, ... % deltau_1min deltau_2min   deltau_1max deltau_2max
                   'RateMinECR',{0,0},'RateMaxECR',{0,0}, ... % V^deltau_1min V^deltau_2min   V^deltau_1max V^deltau_2max (0 = hard constraints)
                   'MinECR',{0,0},'MaxECR',{0,0}, ... % V^u_1min V^u_2min   V^u_1max V^u_2max (0 = hard constraints)
                   'ScaleFactor',{1,1}); % range of u_1 range of u_2     

% Constraints and scaling for OVs (OV = OutputVariables = MO + UO - here we only have MOs)
mpcobj.OV = struct('Min',{10,10},'Max',{90,90}, ... %y1min y2min   y1max y2max
                   'MinECR',{1,1},'MaxECR',{1,1}, ...% V^y_1min V^y_2min   V^y_1max V^y_2max (1 = soft constraints)
                   'ScaleFactor',{80,80}); ... % range of y_1 range of y_2      
% Weights
mpcobj.Weights.OutputVariables = [2 0.5]; % phi1 and phi2
mpcobj.Weights.ManipulatedVariablesRate = [0.1 0.1]; % lambda1 and lambda2
mpcobj.Weights.ECR = 1e5; % rho_epsilon
%

% This should already be set by default
%setoutdist(mpcobj,'integrators');

% Calculate estimator parameters
[Q, R, Gpred] = get_MPC_KF_params_OD(mpcobj);

assert(isequal(round(Q, 6), [ ...
  0  0  0         0         0         0         0  0  0  0  0
  0  0  0         0         0         0         0  0  0  0  0
  0  0  0         0         0         0         0  0  0  0  0
  0  0  0  0.230144  0.223521         0         0  0  0  0  0
  0  0  0  0.223521  0.217088         0         0  0  0  0  0
  0  0  0         0         0  0.848443  0.759178  0  0  0  0
  0  0  0         0         0  0.759178  0.679305  0  0  0  0
  0  0  0         0         0         0         0  0  0  0  0
  0  0  0         0         0         0         0  0  0  0  0
  0  0  0         0         0         0         0  0  0  4  0
  0  0  0         0         0         0         0  0  0  0  4]))
assert(isequal(R, [ ...
     80^2  0
     0     80^2
]))
assert(isequal(round(Gpred.A, 6), [ ...
  0  0  0  0.483333  0  0  0  1  0  0  0
  0  0  1  0         0  0  0  0  0  0  0
  0  0  0  0    0  0  0.6960  0  1  0  0
  0  0  0  0.920044  0  0  0  0  0  0  0
  0  0  0  0  0.866878  0  0  0  0  0  0
  0  0  0  0  0  0.846482  0  0  0  0  0
  0  0  0  0  0  0  0.670320  0  0  0  0
  0  0  0  0  0  0.483333  0  0  0  0  0
  0  0  0  0 -0.497143  0  0  0  0  0  0
  0  0  0  0         0  0  0  0  0  1  0
  0  0  0  0         0  0  0  0  0  0  1
]))
assert(isequal(round(Gpred.B, 6), [ ...
         0         0
         0         0
         0         0
  0.479734         0
  0.465927         0
         0  0.921110
         0  0.824200
         0         0
         0         0
         0         0
         0         0
]))
assert(isequal(round(Gpred.C, 6), [ ...
  1  0  0  0  0  0  0  0  0  80 0
  0  1  0  0  0  0  0  0  0  0  80
]))


%% Test get_MPC_KF_params_ID on SISO example

% Same system from example from documentation for output disturbances
% https://www.mathworks.com/help/mpc/ug/design-estimator-equivalent-to-mpc-built-in-kf.html
x0 = [1.6931; 4.3863];
u0 = 7.7738;
y0 = 0.5;
G = ss([-2 1;-9.7726 -1.3863],[0;1],[0.5 0],0);
Ts = 0.1;
Gd = c2d(G,Ts);

% Inset delays within the discrete state-space matrices
Gd = absorbDelay(Gd); 

% Add an input disturbance to each input:
%      dx/dt = Ax + [B B]*[u; du]
%          y = Cx + [D D]*[u; du]
Gd = ss(Gd.A,[Gd.B Gd.B],Gd.C,[Gd.D Gd.D],Gd.Ts);
Gd = setmpcsignals(Gd,'MV',1,'UD',2,'MO',1);

% Turn off information display
old_status = mpcverbosity('off');

% Define default MPC
mpcobj = mpc(Gd,Ts);

% Set disturbance model to default input disturbance model
% (integrated white noises on each input)
setindist(mpcobj,'integrators'); 

% Set some parameters
mpcobj.Model.Nominal.Y = y0;
mpcobj.Model.Nominal.U = [u0 0]';
mpcobj.Model.Nominal.X = x0;
mpcobj.Weights.OutputVariables = 0.2;

% Calculate estimator parameters
[Q, R, Gpred] = get_MPC_KF_params_ID(mpcobj);
assert(isequal(round(Q, 6), [ ...
    0.000020    0.000408       0
    0.000408    0.008453       0
         0         0        0.01
]))
assert(R == 1)
assert(isequal(round(Gpred.A, 6), [ ...
    0.778227    0.083069  0.004435
   -0.811801    0.829206  0.091939
           0           0         1
]))
assert(isequal(round(Gpred.B, 6), [ ...
    0.004435
    0.091939
           0
]))
assert(isequal(round(Gpred.C, 6), [ ...
    0.5  0   0
]))

% Repeat with MPC which has a continuous-time model
G = ss(G.A,[G.B G.B],G.C,[G.D G.D]);
G = setmpcsignals(G,'MV',1,'UD',2,'MO',1);
mpcobj = mpc(G,Ts);
setindist(mpcobj,'integrators'); 
mpcobj.Model.Nominal.Y = y0;
mpcobj.Model.Nominal.U = [u0 0];
mpcobj.Model.Nominal.X = x0;
mpcobj.Weights.OutputVariables = 0.2;

% Calculate estimator parameters
[Q, R, Gpred] = get_MPC_KF_params_ID(mpcobj);
assert(isequal(round(Q, 6), [ ...
    0.000020    0.000408       0
    0.000408    0.008453       0
         0         0        0.01
]))
assert(R == 1)
assert(isequal(round(Gpred.A, 6), [ ...
    0.778227    0.083069  0.004435
   -0.811801    0.829206  0.091939
           0           0         1
]))
assert(isequal(round(Gpred.B, 6), [ ...
    0.004435
    0.091939
           0
]))
assert(isequal(round(Gpred.C, 6), [ ...
    0.5  0   0
]))

% Return information display settings to previous value
mpcverbosity(old_status);


%% Test get_MPC_KF_params_ID on 2x2 system

% Example system used in MPC_Matlab_ID.mlx
Ts = 0.8; % Sampling period
s = tf('s');
G11 = 0.02*exp(-2*s)/(s*(1+4*s));
G12 = 1.9*exp(-4*s)/(1+6*s);
G21 = -0.74*exp(-6*s)/(1+7*s);
G22 = 0.74*exp(-4*s)/(1+5*s);
G = ss([G11 G12; G21 G22]);
% An input disturbance on each input:
%      dx/dt = Ax + [B B]*[u; du]
%          y = Cx + [D D]*[u; du]
Gcont = ss(G.A,[G.B G.B],G.C,[G.D G.D]);
Gd = c2d(Gcont,Ts,'zoh');
% Inset delays within the discrete state-space matrices
Gd = absorbDelay(Gd);
Gd = minreal(Gd);

mpcverbosity off;
Gd = setmpcsignals(Gd,'MV',1:2,'UD',3:4,'MO',1:2);
mpcobj = mpc(Gd,Ts);

% Setup MPC
% Operating points
u_op = [10;20;0;0];
y_op = [30;40];
mpcobj.model.Nominal.U=u_op;
mpcobj.model.Nominal.Y=y_op;
 
% Prediction and control horizons
mpcobj.PredictionHorizon = 20;
mpcobj.ControlHorizon = 1;

% TODO: Need to get it working with scaling factors on MVs. 'ScaleFactor',{50,50}

% Constraints and scaling for Mvs
mpcobj.MV = struct('Min',{0,0},'Max',{50,50}, ... % u_1min u_2min   u1max u2max
                   'RateMin',{-10,-10},'RateMax',{10,10}, ... % deltau_1min deltau_2min   deltau_1max deltau_2max
                   'RateMinECR',{0,0},'RateMaxECR',{0,0}, ... % V^deltau_1min V^deltau_2min   V^deltau_1max V^deltau_2max (0 = hard constraints)
                   'MinECR',{0,0},'MaxECR',{0,0}, ... % V^u_1min V^u_2min   V^u_1max V^u_2max (0 = hard constraints)
                   'ScaleFactor',{1,1}); % range of u_1 range of u_2     

% Constraints and scaling for OVs (OV = OutputVariables = MO + UO - here we only have MOs)
mpcobj.OV = struct('Min',{10,10},'Max',{90,90}, ... %y1min y2min   y1max y2max
                   'MinECR',{1,1},'MaxECR',{1,1}, ...% V^y_1min V^y_2min   V^y_1max V^y_2max (1 = soft constraints)
                   'ScaleFactor',{80,80}); ... % range of y_1 range of y_2      

% Weights
mpcobj.Weights.OutputVariables = [0.5 0.2]; % phi1 and phi2
mpcobj.Weights.ManipulatedVariablesRate = [0.1 0.1]; % lambda1 and lambda2
mpcobj.Weights.ECR = 1e5; % rho_epsilon

% Set the default input disturbance model
setindist(mpcobj,'integrators');

% Calculate estimator parameters
[Q, R, Gpred] = get_MPC_KF_params_ID(mpcobj);
assert(isequal(round(Q, 4), [ ...
    0.0082    0.0017    0.0171         0         0         0         0
    0.0017    0.0004    0.0035         0         0         0         0
    0.0171    0.0035    0.0357         0         0         0         0
         0         0         0    0.1402    0.1384         0         0
         0         0         0    0.1384    0.1366         0         0
         0         0         0         0         0    0.6400         0
         0         0         0         0         0         0    0.6400
]))
assert(isequal(R, [ ...
    6400           0
       0        6400 ...
]))
assert(isequal(round(Gpred.A, 4), [ ...

    0.8187         0         0         0         0    0.0906         0
    0.3625    1.0000         0         0         0    0.0187         0
         0         0    0.8920         0         0    0.1890         0
         0         0         0    0.8752         0         0    0.3745
         0         0         0         0    0.8521         0    0.3696
         0         0         0         0         0    1.0000         0
         0         0         0         0         0         0    1.0000
]))
assert(isequal(round(Gpred.B, 4), [ ...
    0.0906         0
    0.0187         0
    0.1890         0
         0    0.3745
         0    0.3696
         0         0
         0         0
]))
assert(isequal(round(Gpred.C, 4), [ ...
         0    0.0800         0    0.6333         0         0         0
         0         0   -0.4229         0    0.2960         0         0
]))



%% Test MPC design using Q and R

% Example system used in MPC_Matlab_ID.mlx
Ts = 0.8; % Sampling period
s = tf('s');
G11 = 0.02*exp(-2*s)/(s*(1+4*s));
G12 = 1.9*exp(-4*s)/(1+6*s);
G21 = -0.74*exp(-6*s)/(1+7*s);
G22 = 0.74*exp(-4*s)/(1+5*s);
G = ss([G11 G12; G21 G22]);
% An input disturbance on each input:
%      dx/dt = Ax + [B B]*[u; du]
%          y = Cx + [D D]*[u; du]
Gcont = ss(G.A,[G.B G.B],G.C,[G.D G.D]);
model = c2d(Gcont,Ts,'zoh');
% Inset delays within the discrete state-space matrices
model = absorbDelay(model);
model = minreal(model);

% Estimator model dimensions
n = length(model.A);
nu = size(model.B,2)/2; % number of MVs
ny = size(model.C,1);

mpcverbosity off;
model = setmpcsignals(model,'MV',1:2,'UD',3:4,'MO',1:2);
mpcobj = mpc(model,Ts);

% Augmented model with integrators
Aaug = [model.A model.B(:,1:nu); zeros(nu,n) eye(nu)];
Baug = [model.B(:,1:nu); zeros(nu,nu)];
Caug = [model.C zeros(ny,nu)];
Daug = [zeros(ny,nu)];

% Matrices G and H: see help kalman 
G = eye(size(Baug,1));
H = zeros(size(Caug,1),size(Baug,1));

% Choose covariance matrices
% Measurement noise
R = 0.1*eye(ny);
% Process noise
Q = blkdiag(0.1*eye(n-nu),0.001,0.001);      

% Calculate correction gain (prediction form)
[~,Kpred] = kalman(ss(Aaug,[Baug G],Caug,[Daug H],Ts),Q,R);
setEstimator(mpcobj,Kpred);

%TODO: This example doesn't work. L1 ~= L
%[Q1, R1, Gpred] = get_MPC_KF_params_ID(mpcobj);
% assert(isequal(round(Q1, 4), [ ...
%     0.0082    0.0017    0.0171         0         0         0         0
%     0.0017    0.0004    0.0035         0         0         0         0
%     0.0171    0.0035    0.0357         0         0         0         0
%          0         0         0    0.1402    0.1384         0         0
%          0         0         0    0.1384    0.1366         0         0
%          0         0         0         0         0    0.6400         0
%          0         0         0         0         0         0    0.6400
% ]))
% assert(isequal(R1, [ ...
%     6400           0
%        0        6400 ...
% ]))


%% Test MPC design with odd-shaped system

rng(0)
n = 6;
ny = 2;
nu = 3;
Ts = 0.5;
Gd = drss(n,ny,nu);
Gd.Ts = Ts;
Gd.D(:, :) = 0;  % no direct transmission
Gd = absorbDelay(Gd);

mpcverbosity off;
model = setmpcsignals(Gd,'MV',1:nu,'MO',1:ny);
mpcobj = mpc(model,Ts);
u_op = [10;20;30];
y_op = [5;15];
mpcobj.model.Nominal.U = u_op;
mpcobj.model.Nominal.Y = y_op;

% Custom measurement noise model
mpcobj.Model.Noise = ss(diag([5 10]),'Ts',Ts);

setEstimator(mpcobj,'default')

[Q, R, Gpred] = get_MPC_KF_params_OD(mpcobj);
assert(isequal(round(Q, 4), [ ...
    1.8155   -1.4890    0.6373   -0.0460    0.3827   -0.0529         0         0
   -1.4890    3.4479    0.9008    1.1289   -4.4347    0.9182         0         0
    0.6373    0.9008    1.3417    0.1263   -2.9550    0.4570         0         0
   -0.0460    1.1289    0.1263    2.0172   -0.8146    0.6535         0         0
    0.3827   -4.4347   -2.9550   -0.8146    8.7024   -1.4469         0         0
   -0.0529    0.9182    0.4570    0.6535   -1.4469    0.3790         0         0
         0         0         0         0         0         0    0.2500         0
         0         0         0         0         0         0         0    0.2500
]))
assert(isequal(R, [ ...
    5^2  0
    0    10^2
]))


%% Test MPC design with odd-shaped system

rng(0)
n = 6;
ny = 2;
nu = 3;
nw = 2;
Ts = 0.5;
Gd = drss(n,ny,nu+nw);
Gd.Ts = Ts;
Gd.D(:, :) = 0;  % no direct transmission
Gd = absorbDelay(Gd);

mpcverbosity off;
model = setmpcsignals(Gd,'MV',1:nu,'UD',nu+1:nu+nw,'MO',1:ny);
mpcobj = mpc(model,Ts);
u_op = [10;20;30];
y_op = [5;15];
mpcobj.model.Nominal.U = [u_op; zeros(nw,1)];
mpcobj.model.Nominal.Y = y_op;

% Custom measurement noise model
mpcobj.Model.Noise = ss(diag([5 10]),'Ts',Ts);

setindist(mpcobj,'integrators');

[Q, R, Gpred] = get_MPC_KF_params_ID(mpcobj);
assert(isequal(round(Q, 4), [ ...
    1.2559   -0.1282   -0.0274         0    1.4318         0         0         0
   -0.1282    2.2248    1.1072         0   -3.5059         0         0         0
   -0.0274    1.1072    1.3417   -0.6797   -2.9550         0         0         0
         0         0   -0.6797    0.5850    1.0725         0         0         0
    1.4318   -3.5059   -2.9550    1.0725    8.7024         0         0         0
         0         0         0         0         0         0         0         0
         0         0         0         0         0         0    0.2500         0
         0         0         0         0         0         0         0    0.2500
]))
assert(isequal(R, [ ...
    5^2  0
    0    10^2
]))
