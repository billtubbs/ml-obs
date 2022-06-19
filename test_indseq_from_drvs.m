% Test function indseq_from_drvs.m

assert(indseq_from_drvs(0) == 0)
assert(indseq_from_drvs(false) == 0)

assert(indseq_from_drvs(1) == 1)
assert(indseq_from_drvs(true) == 1)

assert(isequal(indseq_from_drvs([1 0 0]'), [1 0 0]))
assert(isequal(indseq_from_drvs([true false false]'), [1 0 0]))

% Combinations of 2 binary variables
G = [0     1
     0     1
     1     0
     1     1
     0     0];
assert(isequal(indseq_from_drvs(G), [2 2 1 3 0]))

% All combinations of 3 binary values
G = [0 0 0
     1 0 0
     0 1 0
     1 1 0
     0 0 1
     1 0 1
     0 1 1
     1 1 1];
assert(isequal(indseq_from_drvs(G), 0:7))
