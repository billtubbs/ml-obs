function p = prob_gamma_given_w(gamma_k, wk, epsilon, sigma_w)
    % Calcualte using Bayes rule
    p = prob_w_given_gamma(wk, gamma_k, sigma_w) ...
        * prob_gamma(gamma_k, epsilon) ./ prob_w(wk, epsilon, sigma_w);
return