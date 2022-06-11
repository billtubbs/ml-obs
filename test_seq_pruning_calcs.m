% Test pruning method

clear all; clc
rng(2)

% Worked example

% Initialize variables
obj.n_filt = 7;
obj.n_main = 3;
obj.n_hold = 4;
obj.f_main = int16(1:obj.n_main);
obj.f_hold = int16(obj.n_main+1:obj.n_main+obj.n_hold);
obj.p_seq_g_Yk = [0.5 0.15 0.2 0.1 0 0.05 0];
nw = 2;

% Random shuffle
obj.f_main = obj.f_main(randperm(length(obj.f_main)));
obj.f_hold = obj.f_hold(randperm(length(obj.f_hold)));

fprintf("[%d %d %d %d][%d %d %d]\n", obj.f_hold, obj.f_main)

p_seq_g_Yk = obj.p_seq_g_Yk

% Index of current most likely sequence
[~, f_max] = max(obj.p_seq_g_Yk)

% Consistency checks - can be removed later
assert(size(obj.f_hold, 2) == obj.n_hold);
assert(size(obj.f_main, 2) == obj.n_main);
assert(isequal(sort(unique([obj.f_main obj.f_hold])), 1:obj.n_filt));

% Right-shift all filters in holding group. This causes
% the last nw values to 'roll-over' to the left of f_hold.
% e.g. circshift([1 2 3], 1) -> [3 1 2]
obj.f_hold = circshift(obj.f_hold, nw);

disp("After shifting holding group:")
fprintf("[%d %d %d %d][%d %d %d]\n", obj.f_hold, obj.f_main)

% Filters to be moved out of holding group
f_move = obj.f_hold(1:nw)

% Rank hypotheses in main group according to 
% conditional probabilities
[~, i_rank] = sort(obj.p_seq_g_Yk(obj.f_main))

% Select those with lowest probability for pruning
f_to_prune = obj.f_main(i_rank(1:nw))
obj.f_main(i_rank(1:nw)) = f_move;

disp("After moving to main group:")
fprintf("[%d %d %d %d][%d %d %d]\n", obj.f_hold, obj.f_main)

% Clone best sequence and put in holding group:
fprintf("Clone of filters{%d} -> filters{%d}, filters{%d}\n", f_max, f_to_prune)
obj.f_hold(1:nw) = f_to_prune;

disp("After replacing holding group with clone of best:")
fprintf("[%d %d %d %d][%d %d %d]\n", obj.f_hold, obj.f_main)

