function [n, nu, ny] = check_dimensions(A, B, C, D)
% [n, nu, ny] = check_dimensions(A, B, C, D)
% Checks and returns dimensions of state-space system
% represented by matrices A, B, C, D.
%
    n = size(A, 1);
    assert(size(A, 2) == n)
    nu = size(B, 2);
    assert(size(B, 1) == n)
    ny = size(C, 1);
    assert(size(C, 2) == n)
    assert(isequal(size(D), [ny nu]))
end