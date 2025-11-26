function LyapunovInfo = Lyapunov_info(sysInfo, potentialInfo, lambda)


d       = sysInfo.d;
gamma   = sysInfo.gamma;
u       = potentialInfo.u;
R       = potentialInfo.R;  % You may define R explicitly if not present


Id = eye(d);

A = (1 / gamma^2) * u * R + 0.5 * (1 - 2 * lambda)^2 * Id;
B = (1 - 2 * lambda) / gamma * Id;
C = (1 / gamma^2) * Id;

LyapunovInfo.A          = A;
LyapunovInfo.B          = B;
LyapunovInfo.C          = C;
LyapunovInfo.lambda     = lambda;
end
