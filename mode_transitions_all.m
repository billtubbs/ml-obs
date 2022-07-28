function [gamma_km1, gamma_k] = mode_transitions_all(nj)
% Returns two vectors representing all possible mode
% transitions between time k-1 and k for a system
% with nj modes.
%
% Example:
% >> [gamma_km1, gamma_k] = mode_transitions_all(2)
% 
% gamma_km1 =
% 
%      0
%      1
%      0
%      1
% 
% 
% gamma_k =
% 
%      0
%      0
%      1
%      1
%

    [gamma_k, gamma_km1] = meshgrid(0:nj-1, 0:nj-1);
    gamma_k = reshape(gamma_k, [], 1);
    gamma_km1 = reshape(gamma_km1, [], 1);

end