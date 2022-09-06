function [X, Y, Ym] = run_simulation_sys(models,U,V,Gamma,nT)
% Simulate switching system

    [n, ~, ny] = check_dimensions(models{1}.A, models{1}.B, models{1}.C);
    X = zeros(nT+1,n);
    Y = zeros(nT+1,ny);
    Ym = zeros(nT+1,ny);
    xk = zeros(n,1);

    for i = 1:nT+1

        % Switch system
        j = Gamma(i) + 1;

        % Inputs
        uk = U(i,:)';

        % Compute y(k)
        yk = models{j}.C * xk;
        yk_m = yk + V(i);

        % Store results
        X(i, :) = xk';
        Y(i, :) = yk';
        Ym(i, :) = yk_m';

        % Compute x(k+1)
        xk = models{j}.A * xk + models{j}.B * uk;

    end
end