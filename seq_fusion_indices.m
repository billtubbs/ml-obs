function [idx_branch, idx_modes, idx_merge] = seq_fusion_indices(seq, nj)

    % Number of merged hypotheses modelled
    nh = size(seq, 1);

    % Construct the indices for branching and prediction steps
    % for given sequence model
    idx_branch = reshape(repmat((1:nh), 2, 1), [], 1);
    idx_modes = reshape(repmat((0:nj-1)', 1, nh), [], 1);
    seq_kmf_to_k = [seq(idx_branch, :) idx_modes];
    
    % Drop first (oldest) values in sequence
    seq_kmfp1_to_k = seq_kmf_to_k(:, 2:end);
    
    % Drop invalid sequences (not found in seq)
    seq_to_keep = ismember(seq_kmfp1_to_k, seq, 'rows');
    seq_kmfp1_to_k = seq_kmfp1_to_k(seq_to_keep, :);
    idx_branch = idx_branch(seq_to_keep, :);
    idx_modes = idx_modes(seq_to_keep, :);
    
    % Index of identical sequences (rows of seq_kmfp1_to_k)
    % which should be merged
    [matches, idx_merge] = ismember(seq_kmfp1_to_k, seq, 'rows');
    
    assert(all(matches))

end