function seq = indseq_from_drvs(Gamma)
    nd = size(Gamma, 2);
    seq = sum(Gamma .* 2.^(0:nd-1), 2)';
end