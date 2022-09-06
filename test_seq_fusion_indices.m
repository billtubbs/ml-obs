% Test the following functions used by the sequence fusion
% multi-model observers:
% 
%  - combinations.m
%  - combinations_lte.m
%  - prob_combs.m
%  - prob_Gammas.m
%  - seq_fusion_indices.m
%  - seq_fusion_indices_inv.m
%  - calc_merge_idx_seq.m
%

clear all


%% Test combination calculation functions
% and prob_Gammas.m

epsilon = [0.99; 0.01];

S = combinations(5, 0);
assert(isequal(cell2mat(S), zeros(1, 5)))
 
S = combinations(5, 1);
assert(isequal(cell2mat(S), eye(5)))

S = combinations(5, 5);
assert(isequal(cell2mat(S), ones(1, 5)))

S = combinations(3, 2);
S_test = ...
     [1     1     0;
      1     0     1;
      0     1     1];
assert(isequal(cell2mat(S), S_test))

S = combinations_lte(3, 3);
S_test = ...
     [0     0     0;
      1     0     0;
      0     1     0;
      0     0     1;
      1     1     0;
      1     0     1;
      0     1     1;
      1     1     1];
assert(isequal(cell2mat(S), S_test))

p = prob_Gammas(S, epsilon);
assert(size(p, 1) == size(S, 1))
assert(abs(p(1) - epsilon(1)^3) < 1e-10)
assert(all(abs(p(2:4) - epsilon(1)^2*epsilon(2)) < 1e-10))
assert(all(abs(p(5:7) - epsilon(1)*epsilon(2)^2) < 1e-10))
assert(abs(p(8) - epsilon(2)^3) < 1e-10)
assert(abs(sum(p) - 1) < 1e-15)

S = combinations_lte(3, 2);
assert(isequal(cell2mat(S), S_test(1:end-1,:)))

S = combinations_lte(5, 2);
p = prob_Gammas(S, epsilon);
assert(abs(1 - sum(p)) > 1e-10)

assert(isequal(cell2mat(combinations_lte(5, 0)), zeros(1, 5)))

S = combinations_gte_prob(3, 0, epsilon);
assert(isequal(cell2mat(S), zeros(1, 3)))

S = combinations_gte_prob(3, 0.97, epsilon);
assert(isequal(cell2mat(S), zeros(1, 3)))

S = combinations_gte_prob(3, 0.98, epsilon);
assert(isequal(cell2mat(S), [zeros(1, 3); eye(3)]))

S1 = combinations_gte_prob(5, 1, epsilon);
S2 = combinations_lte(5, 5);
assert(isequal(cell2mat(S1), cell2mat(S2)))

% Multi-variable sequences
S = combinations(6, 2);
S = cellfun(@(x) reshape(x, 2, []), S, 'UniformOutput', false);
epsilon = [0.99 0.95; 0.01 0.05];
p = prob_Gammas(S, epsilon);
assert(all(p(1) - 0.01*0.99*0.99*0.05*0.95*0.95 < 1e-10));

% Check for all S
for i = 1:numel(S)
    n_shocks = sum(S{i}, 2);
    n_no_shock = size(S{i},2) - n_shocks;
    prob_calc = prod(epsilon(2,:) .^ (n_shocks')) * ...
        prod(epsilon(1,:) .^ (n_no_shock'));
    assert(abs(prob_calc - p(i)) < 1e-12)
end

f = 3;
n_dist = 2;
m = 1;
d = 5;
S = combinations_lte(f*n_dist, m);
epsilon = [0.01; 0.01];
p = (ones(size(epsilon)) - (ones(size(epsilon)) - epsilon).^d)';
p = [ones(size(p))-p; p];
S = cellfun(@(x) reshape(x, n_dist, []), S, 'UniformOutput', false);
p_seq = prob_Gammas(S, p);
assert(abs(p_seq(1) - (p(1,1)^3*p(1,2)^3)) < 1e-10)
assert(all(abs(p_seq([2 4 6]) - (p(1,1)^2*p(2,1)*p(1,2)^3)) < 1e-10))
assert(all(abs(p_seq([3 5 7]) - (p(1,1)^3*p(1,2)^2*p(2,2))) < 1e-10))


%% Test most_probable_combinations.m

p = [0.95; 0.04; 0.01];
S = most_probable_sequences(p, 1, 0.99);
assert(isequal(S, [1; 2]))

p = [0.01 0.04 0.95];
S = most_probable_sequences(p, 1, 0.99);
assert(isequal(S, [3; 2]))

S = most_probable_sequences(p, 2, 0.99);
S_test = [   3     3;
             3     2;
             2     3;
             3     1];
for i=1:size(S,1)
    assert(isempty(setdiff(S(i,:), S_test(i,:))))
end

S = most_probable_sequences(p, 3, 0.95);
S_test = [   3     3     3;
             3     3     2;
             3     2     3];
for i=1:size(S,1)
    assert(isempty(setdiff(S(i,:), S_test(i,:))))
end


%% Test calc_merge_idx_seq.m

% Test sequence 1
seq = {[0 0 0];
       [1 0 0];
       [0 1 0];
       [0 0 1]};
merge_idx_seq = calc_merge_idx_seq(seq);

test_merge_idx_seq = ...
    [1   1   1
     1   2   2
     2   1   3
     3   3   1];
assert(isequal(merge_idx_seq, test_merge_idx_seq))

% Test sequence 2
seq = {[0 0 0];
       [1 0 0];
       [0 1 0];
       [0 0 1];
       [1 1 0];
       [1 0 1];
       [0 1 1]};
merge_idx_seq = calc_merge_idx_seq(seq);

test_merge_idx_seq = ...
    [1   1   1;
     1   2   2;
     2   1   3;
     3   3   1;
     2   2   4;
     3   4   2;
     4   3   3];
assert(isequal(merge_idx_seq, test_merge_idx_seq))

% Test sequence 3
seq = {[0 0 0 0 0];
       [1 0 0 0 0];
       [0 1 0 0 0];
       [0 0 1 0 0];
       [0 0 0 1 0];
       [0 0 0 0 1];
       [1 1 0 0 0];
       [1 0 1 0 0];
       [1 0 0 1 0];
       [1 0 0 0 1];
       [0 1 1 0 0];
       [0 1 0 1 0];
       [0 1 0 0 1];
       [0 0 1 1 0];
       [0 0 1 0 1];
       [0 0 0 1 1]};
merge_idx_seq = calc_merge_idx_seq(seq);

test_merge_idx_seq = ...
   [1    1    1    1    1;
    1    2    2    2    2;
    2    1    3    3    3;
    3    3    1    4    4;
    4    4    4    1    5;
    5    5    5    5    1;
    2    2    6    6    6;
    3    6    2    7    7;
    4    7    7    2    8;
    5    8    8    8    2;
    6    3    3    9    9;
    7    4    9    3   10;
    8    5   10   10    3;
    9    9    4    4   11;
   10   10    5   11    4;
   11   11   11    5    5];
assert(isequal(merge_idx_seq, test_merge_idx_seq))

% Test sequence 4
seq = {[0 0 0];
       [2 0 0];
       [1 0 0];
       [0 2 0];
       [0 1 0];
       [0 0 2];
       [0 0 1];
       [3 0 0];
       [2 2 0];
       [2 1 0];
       [2 0 2];
       [2 0 1];
       [1 2 0];
       [1 1 0];
       [1 0 2];
       [1 0 1];
       [0 3 0];
       [0 2 2];
       [0 2 1];
       [0 1 2];
       [0 1 1];
       [0 0 3]};
merge_idx_seq = calc_merge_idx_seq(seq);

test_merge_idx_seq = ...
    [1    1    1;
     1    2    2;
     1    3    3;
     2    1    4;
     3    1    5;
     4    4    1;
     5    5    1;
     1    6    6;
     2    2    7;
     3    2    8;
     4    7    2;
     5    8    2;
     2    3    9;
     3    3   10;
     4    9    3;
     5   10    3;
     6    1   11;
     7    4    4;
     8    5    4;
     9    4    5;
    10    5    5;
    11   11    1];
assert(isequal(merge_idx_seq, test_merge_idx_seq))


%% Example 1

nj = 2;  % Number of system modes
f = 3;  % Fusion horizon
m = 3;  % Maximum number of shocks

% Generate all sequences
seq = cell2mat(combinations_lte(f, m));
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

% Generate branching, transition and merge indices - invariant version
[idx_branch, idx_modes, idx_merge] = seq_fusion_indices_inv(seq, nj);
assert(isequal([idx_branch idx_modes idx_merge], test_result))

% Generate branching, transition and merge indices - general version
[idx_branch, idx_modes, idx_merge] = seq_fusion_indices(seq, nj);
% These should contain repeating copies of above
assert(isequal(idx_branch, repmat({test_result(:, 1)}, 1, f)))
assert(isequal(idx_modes, repmat({test_result(:, 2)}, 1, f)))
assert(isequal(cell2mat(idx_merge), [ ...
     1     1     1
     2     3     4
     1     2     2
     2     5     6
     3     1     3
     5     3     7
     4     4     1
     6     7     4
     3     2     5
     5     5     8
     4     6     2
     6     8     6
     7     4     3
     8     7     7
     7     6     5
     8     8     8
]))


%% Example 2

nj = 2;  % Number of system modes
f = 3;  % Fusion horizon
m = 2;  % Maximum number of shocks

% Generate all sequences
seq = cell2mat(combinations_lte(f, m));
assert(isequal(seq, [ ...
   0   0   0
   1   0   0
   0   1   0
   0   0   1
   1   1   0
   1   0   1
   0   1   1 ...
]))

% Generate branching, transition and merge indices - invariant version
[idx_branch, idx_modes, idx_merge] = seq_fusion_indices_inv(seq, nj);
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

% Generate branching, transition and merge indices - general version
[idx_branch, idx_modes, idx_merge] = seq_fusion_indices(seq, nj);
assert(isequal([idx_branch{1} idx_modes{1} idx_merge{1}], [ ...
     1     0     1
     1     1     2
     2     0     1
     2     1     2
     3     0     3
     3     1     5
     4     0     4
     4     1     6
     5     0     3
     5     1     5
     6     0     4
     6     1     6
     7     0     7
]))
assert(isequal([idx_branch{2} idx_modes{2} idx_merge{2}], [ ...
     1     0     1
     1     1     3
     2     0     2
     2     1     5
     3     0     1
     3     1     3
     4     0     4
     4     1     7
     5     0     2
     5     1     5
     6     0     6
     7     0     4
     7     1     7
]))
assert(isequal([idx_branch{3} idx_modes{3} idx_merge{3}], [ ...
     1     0     1
     1     1     4
     2     0     2
     2     1     6
     3     0     3
     3     1     7
     4     0     1
     4     1     4
     5     0     5
     6     0     2
     6     1     6
     7     0     3
     7     1     7
]))


%% Example 3

nj = 2;  % Number of system modes
f = 3;  % Fusion horizon
m = 1;  % Maximum number of shocks

% Generate all sequences
seq = cell2mat(combinations_lte(f, m));
assert(isequal(seq, [ ...
   0   0   0
   1   0   0
   0   1   0
   0   0   1 ...
]))

test_result = [
     1     0     1
     1     1     4
     2     0     1
     2     1     4
     3     0     2
     4     0     3
];

% Generate branching, transition and merge indices - invariant version
[idx_branch, idx_modes, idx_merge] = seq_fusion_indices_inv(seq, nj);
assert(isequal([idx_branch idx_modes idx_merge], test_result))

% Generate branching, transition and merge indices - general version
[idx_branch, idx_modes, idx_merge] = seq_fusion_indices(seq, nj);
assert(isequal([idx_branch{1} idx_modes{1} idx_merge{1}], [ ...
     1     0     1
     1     1     2
     2     0     1
     2     1     2
     3     0     3
     4     0     4
]))
assert(isequal([idx_branch{2} idx_modes{2} idx_merge{2}], [ ...
     1     0     1
     1     1     3
     2     0     2
     3     0     1
     3     1     3
     4     0     4
]))
assert(isequal([idx_branch{3} idx_modes{3} idx_merge{3}], [ ...
     1     0     1
     1     1     4
     2     0     2
     3     0     3
     4     0     1
     4     1     4
]))


%% Example 4 - with detection intervals (1995 version)

% Parameters of SF observer
Q0 = [0.01    0
         0    0];
Bw = [0  1]';
epsilon = 0.01;
sigma_wp = [0.0100    1.0000];
f = 15;
m = 1;
d = 5;
nw = 1;
assert(rem(f, d) == 0, "detection interval not compatible")
n_di = f / d;
nj = 2;

% Construct process noise covariance matrices and switching
% sequences over the fusion horizon, and the prior 
% probabilities of each sequence.
[Q, p_gamma, S] = construct_Q_model_SF95(Q0, Bw, epsilon, ...
    sigma_wp, n_di, m, nw);

% Expand sequences by inserting zeros between times
% when shocks occur.
n_filt = size(S, 1);
seq = cell(n_filt, 1);
for i = 1:n_filt
    seq{i} = int16(zeros(size(S{i}, 1), f));
    seq{i}(:, 1:d:f) = S{i};
    % Alternatively, at end of each detection interval
    %seq{i}(:, d:d:f) = S{i};
end
assert(isequal(seq, {...
    [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]
    [1 0 0 0 0 0 0 0 0 0 0 0 0 0 0]
    [0 0 0 0 0 1 0 0 0 0 0 0 0 0 0]
    [0 0 0 0 0 0 0 0 0 0 1 0 0 0 0] ...
}))

seq = cell2mat(seq);

% Generate sets of branching, mode transition and merge
% indices for each time-step of the sequence
[idx_branch, idx_modes, idx_merge] = seq_fusion_indices(seq, nj);

% Step 1
assert(isequal([idx_branch{1} idx_modes{1} idx_merge{1}], [ ...
     1     0     1
     1     1     2
     2     0     1
     2     1     2
     3     0     3
     4     0     4 ...
]))

% Step 6
assert(isequal([idx_branch{6} idx_modes{6} idx_merge{6}], [ ...
     1     0     1
     1     1     3
     2     0     2
     3     0     1
     3     1     3
     4     0     4 ...
]))

% Step 11
assert(isequal([idx_branch{11} idx_modes{11} idx_merge{11}], [ ...
     1     0     1
     1     1     4
     2     0     2
     3     0     3
     4     0     1
     4     1     4 ...
]))

% No branching or merging during transitions steps between shocks
for i = [2 3 4 5 7 8 9 10 12 13 14 15]
    assert(isequal([idx_branch{i} idx_modes{i} idx_merge{i}], [ ...
         1     0     1
         2     0     2
         3     0     3
         4     0     4 ...
    ]))
end


%% Example 5 - with detection intervals (1998 version)

% Parameters of SF observer
Q0 = [0.01    0
         0    0];
Bw = [0  1]';
epsilon = 0.01;
sigma_wp = [0.0100    1.0000];
f = 15;
m = 1;
d = 5;
nw = 1;
assert(rem(f, d) == 0, "detection interval not compatible")
n_di = f / d;
nj = 2;

% Construct process noise covariance matrices and switching
% sequences over the fusion horizon, and the prior 
% probabilities of each sequence.
alpha = (1 - (1 - epsilon).^d);
var_wp = sigma_wp.^2;
var_wp(2) = var_wp(2) ./ d;
[Q, p_gamma, S] = construct_Q_model_SF(Q0, Bw, alpha, ...
    var_wp, n_di, m, nw);

% Expand sequences by inserting zeros between times
% when shocks occur.
n_filt = size(S, 1);
seq = cell2mat(S);
seq = cell2mat(arrayfun( ...
    @(i) repmat(seq(:, i), 1, d), 1:n_di, 'UniformOutput', false ...
));
assert(isequal(seq, [...
   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0
   1   1   1   1   1   0   0   0   0   0   0   0   0   0   0
   0   0   0   0   0   1   1   1   1   1   0   0   0   0   0
   0   0   0   0   0   0   0   0   0   0   1   1   1   1   1 ...
]))

% Generate sets of branching, mode transition and merge
% indices for each time-step of the sequence
[idx_branch, idx_modes, idx_merge] = seq_fusion_indices(seq, nj);

% This is not a suitable switching sequence since no branching
% or merging occurs and all steps are identical.  Instead, a
% different approach is needed where GPB algorithm is only
% applied at transitions between the detection intervals.
assert(all(cell2mat(idx_branch) == idx_branch{1}, [1 2]))
assert(all(cell2mat(idx_modes) == seq, [1 2]))
assert(all(cell2mat(idx_merge) == idx_merge{1}, [1 2]))

% To achieve the desired results, the fusion horizon must
% be set to a shorter length
nh = 11;
[idx_branch, idx_modes, idx_merge] = seq_fusion_indices(seq, nj, nh);

% Step 1
assert(isequal([idx_branch{1} idx_modes{1} idx_merge{1}], [ ...
     1     0     1
     1     1     2
     2     0     1
     2     1     2
     3     0     3
     4     0     4 ...
]))

% Step 6
assert(isequal([idx_branch{6} idx_modes{6} idx_merge{6}], [ ...
     1     0     1
     1     1     3
     2     0     2
     3     0     1
     3     1     3
     4     0     4 ...
]))

% Step 11
assert(isequal([idx_branch{11} idx_modes{11} idx_merge{11}], [ ...
     1     0     1
     1     1     4
     2     0     2
     3     0     3
     4     0     1
     4     1     4 ...
]))

% No branching or merging during transitions steps between shocks
% During 1st detection interval
for i = [2 3 4 5]
    assert(isequal([idx_branch{i} idx_modes{i} idx_merge{i}], [ ...
         1     0     1
         2     1     2
         3     0     3
         4     0     4 ...
    ]))
end

% During 2nd detection interval
for i = [7 8 9 10]
    assert(isequal([idx_branch{i} idx_modes{i} idx_merge{i}], [ ...
         1     0     1
         2     0     2
         3     1     3
         4     0     4 ...
    ]))
end

% During 3rd detection interval
for i = [12 13 14 15]
    assert(isequal([idx_branch{i} idx_modes{i} idx_merge{i}], [ ...
         1     0     1
         2     0     2
         3     0     3
         4     1     4 ...
    ]))
end


%% Example 6 - more than 2 modes

% Parameters of SF observer
Q0 = [ ...
    0.0100         0         0         0
         0    0.0100         0         0
         0         0         0         0
         0         0         0         0
];
Bw = [ ...
     0     0
     0     0
     1     0
     0     1
];
alpha = [0.0248, 0.0248]';
var_wp = [
    0.1    0.2
    0.1    0.2
];
f = 3;  % No. of detection intervals in fusion horizon
m = 1;  % Maximum number of shocks during fusion horizon
nw = 2;  % Number of indenendent random shock inputs

% Generate all sequences
[Q, p_gamma, seq] = construct_Q_model_SF(Q0, Bw, alpha, ...
                var_wp, f, m, nw);
seq = cell2mat(seq);

% Number of modes
nj = 3;

assert(isequal(seq, [ ...
     0     0     0
     2     0     0
     1     0     0
     0     2     0
     0     1     0
     0     0     2
     0     0     1
]))

test_result = [
     1     0     1
     1     1     7
     1     2     6
     2     0     1
     2     1     7
     2     2     6
     3     0     1
     3     1     7
     3     2     6
     4     0     2
     5     0     3
     6     0     4
     7     0     5
];

% Generate branching, transition and merge indices - invariant version
[idx_branch, idx_modes, idx_merge] = seq_fusion_indices_inv(seq, nj);
assert(isequal([idx_branch idx_modes idx_merge], test_result))

% Generate branching, transition and merge indices - general version
[idx_branch, idx_modes, idx_merge] = seq_fusion_indices(seq, nj);
assert(isequal([idx_branch{1} idx_modes{1} idx_merge{1}], [ ...
     1     0     1
     1     1     3
     1     2     2
     2     0     1
     2     1     3
     2     2     2
     3     0     1
     3     1     3
     3     2     2
     4     0     4
     5     0     5
     6     0     6
     7     0     7
]))
assert(isequal([idx_branch{2} idx_modes{2} idx_merge{2}], [ ...
     1     0     1
     1     1     5
     1     2     4
     2     0     2
     3     0     3
     4     0     1
     4     1     5
     4     2     4
     5     0     1
     5     1     5
     5     2     4
     6     0     6
     7     0     7
]))
assert(isequal([idx_branch{3} idx_modes{3} idx_merge{3}], [ ...
     1     0     1
     1     1     7
     1     2     6
     2     0     2
     3     0     3
     4     0     4
     5     0     5
     6     0     1
     6     1     7
     6     2     6
     7     0     1
     7     1     7
     7     2     6
]))


% Repeat with m = 2
m = 2;  % Maximum number of shocks during fusion horizon
% Generate all sequences
[Q, p_gamma, seq] = construct_Q_model_SF(Q0, Bw, alpha, ...
                var_wp, f, m, nw);
seq = cell2mat(seq);

% Number of modes
nj = 4;

assert(isequal(seq, [ ...
     0     0     0
     2     0     0
     1     0     0
     0     2     0
     0     1     0
     0     0     2
     0     0     1
     3     0     0
     2     2     0
     2     1     0
     2     0     2
     2     0     1
     1     2     0
     1     1     0
     1     0     2
     1     0     1
     0     3     0
     0     2     2
     0     2     1
     0     1     2
     0     1     1
     0     0     3
]))

% Generate branching, transition and merge indices - invariant version
[idx_branch, idx_modes, idx_merge] = seq_fusion_indices_inv(seq, nj);
test_result = [
     1     0     1
     1     1     7
     1     2     6
     1     3    22
     2     0     1
     2     1     7
     2     2     6
     2     3    22
     3     0     1
     3     1     7
     3     2     6
     3     3    22
     4     0     2
     4     1    12
     4     2    11
     5     0     3
     5     1    16
     5     2    15
     6     0     4
     6     1    19
     6     2    18
     7     0     5
     7     1    21
     7     2    20
     8     0     1
     8     1     7
     8     2     6
     8     3    22
     9     0     2
     9     1    12
     9     2    11
    10     0     3
    10     1    16
    10     2    15
    11     0     4
    11     1    19
    11     2    18
    12     0     5
    12     1    21
    12     2    20
    13     0     2
    13     1    12
    13     2    11
    14     0     3
    14     1    16
    14     2    15
    15     0     4
    15     1    19
    15     2    18
    16     0     5
    16     1    21
    16     2    20
    17     0     8
    18     0     9
    19     0    10
    20     0    13
    21     0    14
    22     0    17
];
assert(isequal([idx_branch idx_modes idx_merge], test_result))

% Generate branching, transition and merge indices - general version
[idx_branch, idx_modes, idx_merge] = seq_fusion_indices(seq, nj);
assert(isequal([idx_branch{1} idx_modes{1} idx_merge{1}], [ ...
     1     0     1
     1     1     3
     1     2     2
     1     3     8
     2     0     1
     2     1     3
     2     2     2
     2     3     8
     3     0     1
     3     1     3
     3     2     2
     3     3     8
     4     0     4
     4     1    13
     4     2     9
     5     0     5
     5     1    14
     5     2    10
     6     0     6
     6     1    15
     6     2    11
     7     0     7
     7     1    16
     7     2    12
     8     0     1
     8     1     3
     8     2     2
     8     3     8
     9     0     4
     9     1    13
     9     2     9
    10     0     5
    10     1    14
    10     2    10
    11     0     6
    11     1    15
    11     2    11
    12     0     7
    12     1    16
    12     2    12
    13     0     4
    13     1    13
    13     2     9
    14     0     5
    14     1    14
    14     2    10
    15     0     6
    15     1    15
    15     2    11
    16     0     7
    16     1    16
    16     2    12
    17     0    17
    18     0    18
    19     0    19
    20     0    20
    21     0    21
    22     0    22
]))
assert(isequal([idx_branch{2} idx_modes{2} idx_merge{2}], [ ...
     1     0     1
     1     1     5
     1     2     4
     1     3    17
     2     0     2
     2     1    10
     2     2     9
     3     0     3
     3     1    14
     3     2    13
     4     0     1
     4     1     5
     4     2     4
     4     3    17
     5     0     1
     5     1     5
     5     2     4
     5     3    17
     6     0     6
     6     1    20
     6     2    18
     7     0     7
     7     1    21
     7     2    19
     8     0     8
     9     0     2
     9     1    10
     9     2     9
    10     0     2
    10     1    10
    10     2     9
    11     0    11
    12     0    12
    13     0     3
    13     1    14
    13     2    13
    14     0     3
    14     1    14
    14     2    13
    15     0    15
    16     0    16
    17     0     1
    17     1     5
    17     2     4
    17     3    17
    18     0     6
    18     1    20
    18     2    18
    19     0     7
    19     1    21
    19     2    19
    20     0     6
    20     1    20
    20     2    18
    21     0     7
    21     1    21
    21     2    19
    22     0    22
]))
assert(isequal([idx_branch{3} idx_modes{3} idx_merge{3}], [ ...
     1     0     1
     1     1     7
     1     2     6
     1     3    22
     2     0     2
     2     1    12
     2     2    11
     3     0     3
     3     1    16
     3     2    15
     4     0     4
     4     1    19
     4     2    18
     5     0     5
     5     1    21
     5     2    20
     6     0     1
     6     1     7
     6     2     6
     6     3    22
     7     0     1
     7     1     7
     7     2     6
     7     3    22
     8     0     8
     9     0     9
    10     0    10
    11     0     2
    11     1    12
    11     2    11
    12     0     2
    12     1    12
    12     2    11
    13     0    13
    14     0    14
    15     0     3
    15     1    16
    15     2    15
    16     0     3
    16     1    16
    16     2    15
    17     0    17
    18     0     4
    18     1    19
    18     2    18
    19     0     4
    19     1    19
    19     2    18
    20     0     5
    20     1    21
    20     2    20
    21     0     5
    21     1    21
    21     2    20
    22     0     1
    22     1     7
    22     2     6
    22     3    22
]))


%% Mock simulations 1

% Run a mock simulation to see what sequences are 
% created when initializing with a set of identical
% estimates.

nj = 2;  % Number of system modes
f = 3;  % Fusion horizon
m = 3;  % Maximum number of shocks

% Generate all sequences
seq = cell2mat(combinations_lte(f, m));
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
nh = size(seq, 1);

% Number of transitions to simulate
nT = 3;

% Generate branching, transition and merge indices - invariant version
[idx_branch1, idx_modes1, idx_merge1] = seq_fusion_indices_inv(seq, nj);
sim_seq1 = mock_sim_inv(nh, idx_branch1, idx_modes1, idx_merge1, nT);

% Generate branching, transition and merge indices - general version
[idx_branch2, idx_modes2, idx_merge2] = seq_fusion_indices(seq, nj);
sim_seq2 = mock_sim(nh, idx_branch2, idx_modes2, idx_merge2, nT);

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
assert(isequal(sim_seq1, test_result))
assert(isequal(sortrows(sim_seq2), sortrows(test_result)))

% Number of transitions to simulate
nT = 4;

% Generate branching, transition and merge indices - invariant version
[idx_branch1, idx_modes1, idx_merge1] = seq_fusion_indices_inv(seq, nj);
sim_seq1 = mock_sim_inv(nh, idx_branch1, idx_modes1, idx_merge1, nT);

% Generate branching, transition and merge indices - general version
[idx_branch2, idx_modes2, idx_merge2] = seq_fusion_indices(seq, nj);
sim_seq2 = mock_sim(nh, idx_branch2, idx_modes2, idx_merge2, nT);

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
assert(isequal(sim_seq1, test_result))
assert(isequal(sortrows(sim_seq2), sortrows(test_result)))


%% Mock simulations 2

nj = 2;  % Number of system modes
f = 3;  % Fusion horizon
m = 1;  % Maximum number of shocks

% Generate all sequences
seq = cell2mat(combinations_lte(f, m));
assert(isequal(seq, [ ...
   0   0   0
   1   0   0
   0   1   0
   0   0   1 ...
]))
nh = size(seq, 1);

% Generate branching, mode transition and merge indices
[idx_branch, idx_modes, idx_merge] = seq_fusion_indices_inv(seq, nj);

% Number of transitions to simulate
nT = 3;

% Generate branching, transition and merge indices - invariant version
[idx_branch1, idx_modes1, idx_merge1] = seq_fusion_indices_inv(seq, nj);
sim_seq1 = mock_sim_inv(nh, idx_branch1, idx_modes1, idx_merge1, nT);

% Generate branching, transition and merge indices - general version
[idx_branch2, idx_modes2, idx_merge2] = seq_fusion_indices(seq, nj);
sim_seq2 = mock_sim(nh, idx_branch2, idx_modes2, idx_merge2, nT);

test_result = [
    "x|0|0|0"
    "x|1|0|0"
    "x|0|1|0"
    "x|0|0|1"
];
assert(isequal(sim_seq1, test_result))
assert(isequal(sortrows(sim_seq2), sortrows(test_result)))

% Number of transitions to simulate
nT = 4;

% Generate branching, transition and merge indices - invariant version
[idx_branch1, idx_modes1, idx_merge1] = seq_fusion_indices_inv(seq, nj);
sim_seq1 = mock_sim_inv(nh, idx_branch1, idx_modes1, idx_merge1, nT);

% Generate branching, transition and merge indices - general version
[idx_branch2, idx_modes2, idx_merge2] = seq_fusion_indices(seq, nj);
sim_seq2 = mock_sim(nh, idx_branch2, idx_modes2, idx_merge2, nT);

test_result = [
    "(x|0|0|0|0+x|1|0|0|0)"
    "x|0|1|0|0"
    "x|0|0|1|0"
    "(x|0|0|0|1+x|1|0|0|1)"
];
assert(isequal(sim_seq1, test_result))
assert(isequal(sortrows(sim_seq2), sortrows(test_result)))

nT = 5;

% Generate branching, transition and merge indices - invariant version
[idx_branch1, idx_modes1, idx_merge1] = seq_fusion_indices_inv(seq, nj);
sim_seq1 = mock_sim_inv(nh, idx_branch1, idx_modes1, idx_merge1, nT);

% Generate branching, transition and merge indices - general version
[idx_branch2, idx_modes2, idx_merge2] = seq_fusion_indices(seq, nj);
sim_seq2 = mock_sim(nh, idx_branch2, idx_modes2, idx_merge2, nT);

test_result = [
    "((x|0|0|0|0+x|1|0|0|0)|0+x|0|1|0|0|0)"
    "x|0|0|1|0|0"
    "(x|0|0|0|1+x|1|0|0|1)|0"
    "((x|0|0|0|0+x|1|0|0|0)|1+x|0|1|0|0|1)"
];
assert(isequal(sim_seq1, test_result))
assert(isequal(sortrows(sim_seq2), sortrows(test_result)))


%% Mock simulations 3

seq = cell2mat({...
    [0 0 0 0 0 0 0 0 0]
    [1 0 0 0 0 0 0 0 0]
    [0 0 0 1 0 0 0 0 0]
    [0 0 0 0 0 0 1 0 0] ...
});
nh = size(seq, 1);
nj = 2;

% Generate branching, transition and merge indices - general version
[idx_branch, idx_modes, idx_merge] = seq_fusion_indices(seq, nj);

% Number of transitions to simulate
nT = 9;
sim_seq = mock_sim(nh, idx_branch, idx_modes, idx_merge, nT);

test_result = [
    "x|0|0|0|0|0|0|0|0|0"
    "x|1|0|0|0|0|0|0|0|0"
    "x|0|0|0|1|0|0|0|0|0"
    "x|0|0|0|0|0|0|1|0|0"
];
assert(isequal(sim_seq, test_result))


% Functions used in the mock sim tests

function sim_seq = mock_sim_inv(nh, idx_branch, idx_modes, idx_merge, nT)
% Produces a string representation of the sequences of
% branching, mode transitions, and hypothesis merging
% steps based on given indices.  This version is for the
% case where  idx_branch, idx_modes, and idx_merge were
% generated with the function seq_fusion_indices_inv.
%

    sim_seq = repmat("x", nh, 1);

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


function sim_seq = mock_sim(nh, idx_branch, idx_modes, idx_merge, nT)
% Produces a string representation of the sequences of
% branching, mode transitions, and hypothesis merging
% steps based on given indices.  This version is for the
% case where  idx_branch, idx_modes, and idx_merge were
% generated with the function seq_fusion_indices.
%

    sim_seq = repmat("x", nh, 1);

    % Length of sequences
    f = size(idx_branch, 2);

    for k = 1:nT

        % Location in sequence
        j = mod(k - 1, f) + 1;

        % Branching
        seq_branched = sim_seq(idx_branch{j}, 1);

        % Prediction step
        seq_branched = seq_branched + "|" + string(idx_modes{j});

        % Merging step
        idx_merged = unique(idx_merge{j});
        nm = length(idx_merged);
        for i = 1:nm
            idx = find(idx_merge{j} == i);
            if all(seq_branched(idx) == seq_branched(idx(1)))
                sim_seq(i) = seq_branched(idx(1));
            else
                sim_seq(i) = sprintf("(%s)", strjoin(seq_branched(idx), "+"));
            end
        end
    end
end
