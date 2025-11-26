function r2 = compute_Lyapunov_distance(x, v, x_ref, v_ref, LyapunovInfo)
%COMPUTE_LYAPUNOV_DISTANCE Compute r^2((x,v), (x',v')) from matrix coefficients
%
% Inputs:
%   x, v     - current state vectors (d × 1)
%   x_ref, v_ref - reference state vectors (d × 1)
%   A, B, C  - matrices of size d × d
%
% Output:
%   r2       - scalar, value of r^2

    A = LyapunovInfo.A;
    B = LyapunovInfo.B;
    C = LyapunovInfo.C;

    [~, N] = size(x);
    r2 = zeros(N, 1);
    
    x_diff = x - x_ref;
    v_diff = v - v_ref;

    for n = 1:N
        dx = x_diff(:,n);
        dv = v_diff(:,n);
        r2(n) = dx' * A * dx + dx' * B * dv + dv' * C * dv;
    end


    r2(r2 <=0) = eps;
end