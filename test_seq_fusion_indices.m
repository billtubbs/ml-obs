% Testing sequence fusion index generation functions

clear all


%% Example 1

nj = 2;  % Number of system modes
n = 3;  % Fusion horizon
m = 3;  % Maximum number of shocks

% Generate all sequences
seq = cell2mat(combinations_lte(n, m));
assert(isequal(seq, [ ...
   0   0   0
   1   0   0
   0   1   0
   0   0   1
   1   1   0
   1   0   1
   0   1   1
   1   1   1 ...
]))

% Generate branching, mode transition and merge indices
[idx_branch, idx_modes, idx_merge] = seq_fusion_indices(seq, nj);

% Show results
test_result = [
     1     0     1
     1     1     4
     2     0     1
     2     1     4
     3     0     2
     3     1     6
     4     0     3
     4     1     7
     5     0     2
     5     1     6
     6     0     3
     6     1     7
     7     0     5
     7     1     8
     8     0     5
     8     1     8
];
assert(isequal([idx_branch idx_modes idx_merge], test_result))


%% Example 2

nj = 2;  % Number of system modes
n = 3;  % Fusion horizon
m = 2;  % Maximum number of shocks

% Generate all sequences
seq = cell2mat(combinations_lte(n, m));
assert(isequal(seq, [ ...
   0   0   0
   1   0   0
   0   1   0
   0   0   1
   1   1   0
   1   0   1
   0   1   1 ...
]))

% Generate branching, mode transition and merge indices
[idx_branch, idx_modes, idx_merge] = seq_fusion_indices(seq, nj);

% Show results
test_result = [
     1     0     1
     1     1     4
     2     0     1
     2     1     4
     3     0     2
     3     1     6
     4     0     3
     4     1     7
     5     0     2
     5     1     6
     6     0     3
     6     1     7
     7     0     5
];
assert(isequal([idx_branch idx_modes idx_merge], test_result))


%% Example 3

nj = 2;  % Number of system modes
n = 3;  % Fusion horizon
m = 1;  % Maximum number of shocks

% Generate all sequences
seq = cell2mat(combinations_lte(n, m));
assert(isequal(seq, [ ...
   0   0   0
   1   0   0
   0   1   0
   0   0   1 ...
]))

% Generate branching, mode transition and merge indices
[idx_branch, idx_modes, idx_merge] = seq_fusion_indices(seq, nj);

% Show results
test_result = [
     1     0     1
     1     1     4
     2     0     1
     2     1     4
     3     0     2
     4     0     3
];
assert(isequal([idx_branch idx_modes idx_merge], test_result))


%% Mock simulations 1

% Run a mock simulation to see what sequences are 
% created when initializing with a set of identical
% estimates.

nj = 2;  % Number of system modes
n = 3;  % Fusion horizon
m = 3;  % Maximum number of shocks

% Generate all sequences
seq = cell2mat(combinations_lte(n, m));
nh = size(seq, 1);

% Generate branching, mode transition and merge indices
[idx_branch, idx_modes, idx_merge] = seq_fusion_indices(seq, nj);

nT = 3;
sim_seq = mock_sim(nh, idx_branch, idx_modes, idx_merge, nT);
test_result = [
    "x|0|0|0"
    "x|1|0|0"
    "x|0|1|0"
    "x|0|0|1"
    "x|1|1|0"
    "x|1|0|1"
    "x|0|1|1"
    "x|1|1|1"
];
assert(isequal(sim_seq, test_result))

nT = 4;
sim_seq = mock_sim(nh, idx_branch, idx_modes, idx_merge, nT);
test_result = [
    "(x|0|0|0|0+x|1|0|0|0)"
    "(x|0|1|0|0+x|1|1|0|0)"
    "(x|0|0|1|0+x|1|0|1|0)"
    "(x|0|0|0|1+x|1|0|0|1)"
    "(x|0|1|1|0+x|1|1|1|0)"
    "(x|0|1|0|1+x|1|1|0|1)"
    "(x|0|0|1|1+x|1|0|1|1)"
    "(x|0|1|1|1+x|1|1|1|1)"
];
assert(isequal(sim_seq, test_result))


%% Mock simulations 2

nj = 2;  % Number of system modes
n = 3;  % Fusion horizon
m = 1;  % Maximum number of shocks

% Generate all sequences
seq = cell2mat(combinations_lte(n, m));
nh = size(seq, 1);

% Generate branching, mode transition and merge indices
[idx_branch, idx_modes, idx_merge] = seq_fusion_indices(seq, nj);

nT = 3;
sim_seq = mock_sim(nh, idx_branch, idx_modes, idx_merge, nT);
test_result = [
    "x|0|0|0"
    "x|1|0|0"
    "x|0|1|0"
    "x|0|0|1"
];
assert(isequal(sim_seq, test_result))

nT = 4;
sim_seq = mock_sim(nh, idx_branch, idx_modes, idx_merge, nT);
test_result = [
    "(x|0|0|0|0+x|1|0|0|0)"
    "x|0|1|0|0"
    "x|0|0|1|0"
    "(x|0|0|0|1+x|1|0|0|1)"
];
assert(isequal(sim_seq, test_result))

nT = 5;
sim_seq = mock_sim(nh, idx_branch, idx_modes, idx_merge, nT);
test_result = [
    "((x|0|0|0|0+x|1|0|0|0)|0+x|0|1|0|0|0)"
    "x|0|0|1|0|0"
    "(x|0|0|0|1+x|1|0|0|1)|0"
    "((x|0|0|0|0+x|1|0|0|0)|1+x|0|1|0|0|1)"
];
assert(isequal(sim_seq, test_result))


function sim_seq = mock_sim(nh, idx_branch, idx_modes, idx_merge, nT)
% Produces a string representation of the sequences of
% branching, mode transitions, and hypothesis merging
% steps based on given indices.
%
    sim_seq = repmat("x",nh,1);
    for k = 1:nT
        
        % Branching
        seq_branched = sim_seq(idx_branch, 1);
    
        % Prediction step
        seq_branched = seq_branched + "|" + string(idx_modes);
    
        % Merging step
        idx_merged = unique(idx_merge);
        nm = length(idx_merged);
        for i = 1:nm
            idx = find(idx_merge == i);
            if all(seq_branched(idx) == seq_branched(idx(1)))
                sim_seq(i) = seq_branched(idx(1));
            else
                sim_seq(i) = sprintf("(%s)", strjoin(seq_branched(idx), "+"));
            end
        end
    end
end
