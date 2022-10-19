% Test get_MPC_KF_params_OD.m and get_MPC_KF_params_ID.m


%% Test get_MPC_KF_params_OD

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


%% Test

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
mpcobj.Weights.OutputVariables = [0.5 0.2];

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
     1     0
     0     1
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
  1  0  0  0  0  0  0  0  0  1  0
  0  1  0  0  0  0  0  0  0  0  1
]))


%% Test get_MPC_KF_params_ID

% Based on example from documentation with output disturbances
% https://www.mathworks.com/help/mpc/ug/design-estimator-equivalent-to-mpc-built-in-kf.html

x0 = [1.6931; 4.3863];
u0 = 7.7738;
y0 = 0.5;
G = ss([-2 1;-9.7726 -1.3863],[0;1],[0.5 0],0);
Ts = 0.1;
Gd = c2d(G,Ts);

% Add an input disturbance to each input:
%      dx/dt = Ax + [B B]*[u; du]
%          y = Cx + [D D]*[u; du]
Gd = ss(Gd.A,[Gd.B Gd.B],Gd.C,[Gd.D Gd.D]);
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

% Return information display settings to previous value
mpcverbosity(old_status);