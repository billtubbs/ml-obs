function kalman_filter_sfunc(block)
%% Simulate a standard Kalman filter in Simulink
%
%

setup(block);

%end kalman_filter_sfunc



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

% Input 1: u(k)
block.InputPort(1).Dimensions        = 1;
block.InputPort(1).DatatypeID  = 0;  % double
block.InputPort(1).Complexity  = 'Real';
block.InputPort(1).DirectFeedthrough = false;
block.InputPort(1).SamplingMode = 'Sample';

% Input 2: y(k)
block.InputPort(2).Dimensions        = 1;
block.InputPort(2).DatatypeID  = 0;  % double
block.InputPort(2).Complexity  = 'Real';
block.InputPort(2).DirectFeedthrough = false;
block.InputPort(2).SamplingMode = 'Sample';

% Output 1: x_est(k+1);
block.OutputPort(1).Dimensions       = 2;
block.OutputPort(1).DatatypeID  = 0; % double
block.OutputPort(1).Complexity  = 'Real';
block.OutputPort(1).SamplingMode = 'Sample';

% Output 2: y_est(k+1)
block.OutputPort(2).Dimensions       = 1;
block.OutputPort(2).DatatypeID  = 0; % double
block.OutputPort(2).Complexity  = 'Real';
block.OutputPort(2).SamplingMode = 'Sample';

% Register parameters (number of parameters in the dialog box)
block.NumDialogPrms     = 1;

% Register sample times
%  [0 offset]            : Continuous sample time
%  [positive_num offset] : Discrete sample time
%
%  [-1, 0]               : Inherited sample time
%  [-2, 0]               : Variable sample time
block.SampleTimes = [-1 0];

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

block.RegBlockMethod('PostPropagationSetup',    @DoPostPropSetup);
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

block.NumDworks = 4;

% Dynamic state estimates: xkp1_est
block.Dwork(1).Name            = 'xkp1_est';
block.Dwork(1).Dimensions      = 2;
block.Dwork(1).DatatypeID      = 0;      % double
block.Dwork(1).Complexity      = 'Real'; % real
block.Dwork(1).UsedAsDiscState = true;

% Dynamic output estimates: ykp1_est
block.Dwork(2).Name            = 'ykp1_est';
block.Dwork(2).Dimensions      = 1;
block.Dwork(2).DatatypeID      = 0;      % double
block.Dwork(2).Complexity      = 'Real'; % real
block.Dwork(2).UsedAsDiscState = true;

% Dynamic estimate covariance matrix: P
block.Dwork(3).Name            = 'P';
block.Dwork(3).Dimensions      = 4;
block.Dwork(3).DatatypeID      = 0;      % double
block.Dwork(3).Complexity      = 'Real'; % real
block.Dwork(3).UsedAsDiscState = true;

% Dynamic gain matrix: K
block.Dwork(4).Name            = 'K';
block.Dwork(4).Dimensions      = 2;
block.Dwork(4).DatatypeID      = 0;      % double
block.Dwork(4).Complexity      = 'Real'; % real
block.Dwork(4).UsedAsDiscState = true;

%end PostPropagationSetup



function InitializeConditions(block)
%%   Functionality    : Called at the start of simulation and if it is 
%                      present in an enabled subsystem configured to reset 
%                      states, it will be called when the enabled subsystem
%                      restarts execution to reset the states.
%   Required         : No
%   C MEX counterpart: mdlInitializeConditions  

x0 = [0; 0];
y0 = 0;
P0 = [ 0.0857 0.0444;
       0.0444 0.0368];
K0 = nan(2, 1);

% Initialize Dwork
block.Dwork(1).Data = reshape(x0, 1, []);
block.Dwork(2).Data = reshape(y0, 1, []);
block.Dwork(3).Data = reshape(P0, 1, []);
block.Dwork(4).Data = reshape(K0, 1, []);

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

% Output y_est(k+1)
block.OutputPort(1).Data = block.Dwork(1).Data;

% Output y_est(k+1)
block.OutputPort(2).Data = block.Dwork(2).Data;

%end Outputs



function Update(block)
%%   Functionality    : Called to update discrete states
%                      during simulation step
%   Required         : No
%   C MEX counterpart: mdlUpdate

% Get observer struct
obs = block.DialogPrm(1).Data;

% Inputs
uk = block.InputPort(1).Data;
yk = block.InputPort(2).Data;

% Variables from memory
xkp1_est = reshape(block.Dwork(1).Data, [2 1]);
P = reshape(block.Dwork(3).Data, [2 2]);
K = reshape(block.Dwork(4).Data, [2 1]);

% Calculate Kalman filter updates
[K, P] = kalman_update(P, obs.A, obs.C, obs.Q, obs.R);

% Update state and output estimates for next timestep
xkp1_est = obs.A * xkp1_est + obs.B * uk + K * (yk - obs.C * xkp1_est);
ykp1_est = obs.C * xkp1_est;

% Save updated variables
block.Dwork(1).Data = reshape(xkp1_est, 1, []);
block.Dwork(2).Data = reshape(ykp1_est, 1, []);
block.Dwork(3).Data = reshape(P, 1, []);
block.Dwork(4).Data = reshape(K, 1, []);

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


