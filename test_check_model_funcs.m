% Test model parameter-checking functions
% 
%  - check_dimensions.m
%  - check_model.m
%  - check_models.m
%

clear all
rng(0)

% Dimensions of test systems
p.n = 4; p.ny = 2; p.nu = 3;

% Sampling period
p.Ts = 0.5;

% Number of systems
nj = 3;

% Create random systems with given dimensions
n_models = nj;
models = cell(1, n_models);
for i = 1:n_models
    models{i} = drss(p.n,p.ny,p.nu);
    models{i}.Ts = p.Ts;
end

model = models{1};
[test.n, test.nu, test.ny] = check_dimensions(model.A, model.B, model.C);
assert(isequal([test.n, test.nu, test.ny], [p.n, p.nu, p.ny]))

[test.n, test.nu, test.ny] = check_dimensions(model.A, model.B, model.C, model.D);
assert(isequal([test.n, test.nu, test.ny], [p.n, p.nu, p.ny]))

[test.n, test.nu, test.ny, test.Ts, direct] = check_model(model);
assert(isequal(test, p))

[test_nj, test.n, test.nu, test.ny, test.Ts] = check_models(models);
assert(test_nj == nj)
assert(isequal(test, p))

% TODO Check ValueErrors occur
model_a = drss(2,1,1);

