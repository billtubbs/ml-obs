% Test check_model.m check_models.m, and check_dimensions

n = 4; ny = 2; nu = 2;

n_models = 3;
models = cell(1, n_models);
for i = 1:n_models
    models{i} = rss(n,ny,nu);
end


