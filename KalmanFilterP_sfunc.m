function KalmanFilterP_sfunc(block)
%% Simulate a Kalman filter in Simulink
%
% This version is for the KalmanFilterP class, which is the
% prediction version that calculates the predictions (i.e. 
% prior estimates) of the states and outputs, x_est(k|k-1) 
% and y_est(k|k-1) at time k.
%

setup(block);

%end KalmanFilterP_sfunc



function setup(block)
%%   Set up the basic characteristics of the S-function block such as:
%      - Input ports
%      - Output ports
%      - Dialog parameters
%      - Options
%   Required         : Yes
%   C MEX counterpart: mdlInitializeSizes

% Register number of ports
block.NumInputPorts  = 2;
block.NumOutputPorts = 2;

% Setup port properties to be inherited or dynamic
block.SetPreCompInpPortInfoToDynamic;
block.SetPreCompOutPortInfoToDynamic;

% Get observer object
obs = block.DialogPrm(1).Data;

% Input 1: u(k)
block.InputPort(1).Dimensions = obs.nu;
block.InputPort(1).DatatypeID = 0;  % double
block.InputPort(1).Complexity = 'Real';
block.InputPort(1).DirectFeedthrough = false;
block.InputPort(1).SamplingMode = 'Sample';

% Input 2: y(k)
block.InputPort(2).Dimensions = obs.ny;
block.InputPort(2).DatatypeID = 0;  % double
block.InputPort(2).Complexity = 'Real';
block.InputPort(2).DirectFeedthrough = false;
block.InputPort(2).SamplingMode = 'Sample';

% Output 1: x_est(k|k-1);
block.OutputPort(1).Dimensions = obs.n;
block.OutputPort(1).DatatypeID = 0; % double
block.OutputPort(1).Complexity = 'Real';
block.OutputPort(1).SamplingMode = 'Sample';

% Output 2: y_est(k|k-1)
block.OutputPort(2).Dimensions = obs.ny;
block.OutputPort(2).DatatypeID = 0; % double
block.OutputPort(2).Complexity = 'Real';
block.OutputPort(2).SamplingMode = 'Sample';

% Register parameters (number of parameters in the dialog box)
block.NumDialogPrms = 1;

% Register sample times
%  [0 offset]            : Continuous sample time
%  [positive_num offset] : Discrete sample time
%
%  [-1, 0]               : Inherited sample time
%  [-2, 0]               : Variable sample time
block.SampleTimes = [obs.Ts 0];

% Specify the block simStateCompliance. The allowed values are:
%    'UnknownSimState', < The default setting; warn and assume DefaultSimState
%    'DefaultSimState', < Same sim state as a built-in block
%    'HasNoSimState',   < No sim state
%    'CustomSimState',  < Has GetSimState and SetSimState methods
%    'DisallowSimState' < Error out when saving or restoring the model sim state
block.SimStateCompliance = 'DefaultSimState';

% -----------------------------------------------------------------
% The MATLAB S-function uses an internal registry for all
% block methods. You should register all relevant methods
% (optional and required) as illustrated below. You may choose
% any suitable name for the methods and implement these methods
% as local functions within the same file. See comments
% provided for each function for more information.
% -----------------------------------------------------------------

block.RegBlockMethod('PostPropagationSetup', @DoPostPropSetup);
block.RegBlockMethod('InitializeConditions', @InitializeConditions);
%block.RegBlockMethod('Start', @Start);
block.RegBlockMethod('Outputs', @Outputs);     % Required
block.RegBlockMethod('Update', @Update);
%block.RegBlockMethod('Derivatives', @Derivatives);
block.RegBlockMethod('Terminate', @Terminate); % Required

%end setup


function DoPostPropSetup(block)
%% PostPropagationSetup:
%   Functionality    : Setup work areas and state variables. Can
%                      also register run-time methods here
%   Required         : No
%   C MEX counterpart: mdlSetWorkWidths

% Get observer object
obs = block.DialogPrm(1).Data;

% Determine system dimensions
%[n, nu, ny] = check_dimensions(obs.A, obs.B, obs.C, obs.D);

block.NumDworks = 3;

% Predictions of states in next time instant: x_est(k+1|k)
block.Dwork(1).Name            = 'xkp1_est';
block.Dwork(1).Dimensions      = obs.n;
block.Dwork(1).DatatypeID      = 0;      % double
block.Dwork(1).Complexity      = 'Real'; % real
block.Dwork(1).UsedAsDiscState = true;

% Predictions of outputs in next time instant: y_est(k+1|k)
block.Dwork(2).Name            = 'ykp1_est';
block.Dwork(2).Dimensions      = obs.ny;
block.Dwork(2).DatatypeID      = 0;      % double
block.Dwork(2).Complexity      = 'Real'; % real
block.Dwork(2).UsedAsDiscState = true;

% State prediction error covariance: P(k+1|k)
block.Dwork(3).Name            = 'Pkp1';
block.Dwork(3).Dimensions      = obs.n*obs.n;
block.Dwork(3).DatatypeID      = 0;      % double
block.Dwork(3).Complexity      = 'Real'; % real
block.Dwork(3).UsedAsDiscState = true;

%end PostPropagationSetup



function InitializeConditions(block)
%%   Functionality    : Called at the start of simulation and if it is 
%                      present in an enabled subsystem configured to reset 
%                      states, it will be called when the enabled subsystem
%                      restarts execution to reset the states.
%   Required         : No
%   C MEX counterpart: mdlInitializeConditions  

% Get observer object
obs = block.DialogPrm(1).Data;

% Initialize Dwork
block.Dwork(1).Data = obs.xkp1_est;
block.Dwork(2).Data = obs.ykp1_est;
block.Dwork(3).Data = obs.Pkp1(:);

%end InitializeConditions



%function Start(block)
%%   Functionality    : Called once at start of model execution. If you
%                      have states that should be initialized once, this 
%                      is the place to do it.
%   Required         : No
%   C MEX counterpart: mdlStart%end Start

%end Start



function Outputs(block)
%%   Functionality    : Called to generate block outputs in
%                      simulation step
%   Required         : Yes
%   C MEX counterpart: mdlOutputs

% Output x_est(k|k-1)
block.OutputPort(1).Data = block.Dwork(1).Data;

% Output y_est(k|k-1)
block.OutputPort(2).Data = block.Dwork(2).Data;

%end Outputs



function Update(block)
%%   Functionality    : Called to update discrete states
%                      during simulation step
%   Required         : No
%   C MEX counterpart: mdlUpdate

% Get observer object
obs = block.DialogPrm(1).Data;

% Inputs
uk = block.InputPort(1).Data;
yk = block.InputPort(2).Data;

% Check size of input vectors
assert(isequal(size(uk), [obs.nu 1]))
assert(isequal(size(yk), [obs.ny 1]))

% Variables from memory
xkp1_est = block.Dwork(1).Data;
ykp1_est = block.Dwork(2).Data;
Pkp1 = reshape(block.Dwork(3).Data, [obs.n obs.n]);

% Get observer object
obs = block.DialogPrm(1).Data;

% Update variables from Simulink memory
obs.xkp1_est = xkp1_est;
obs.ykp1_est = ykp1_est;
obs.Pkp1 = reshape(Pkp1, obs.n, obs.n);

% Kalman filter update
obs.update(yk, uk);

% % Calculate Kalman filter updates
% [K, P] = kalman_update(P, obs.A, obs.C, obs.Q, obs.R);
% 
% % Update state and output estimates for next timestep
% xkp1_est = obs.A * xkp1_est + obs.B * uk + K * (yk - ykp1_est);
% ykp1_est = obs.C * xkp1_est;

% Save updated variables as row vectors
block.Dwork(1).Data = obs.xkp1_est;
block.Dwork(2).Data = obs.ykp1_est;
block.Dwork(3).Data = obs.Pkp1(:);

%end Update



%function Derivatives(block)
%%   Functionality    : Called to update derivatives of
%                      continuous states during simulation step
%   Required         : No
%   C MEX counterpart: mdlDerivatives

%end Derivatives



function Terminate(block)
%%   Functionality    : Called at the end of simulation for cleanup
%   Required         : Yes
%   C MEX counterpart: mdlTerminate

%end Terminate
