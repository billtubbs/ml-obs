function p = prob_w_given_gamma(wk, gamma_k, sigma_w)
    switch gamma_k
        case 0
            p = normpdf(wk, 0, sigma_w(1));
        case 1
            p = normpdf(wk, 0, sigma_w(2));
        otherwise
            error("Value error: gamma_k")
    end
return