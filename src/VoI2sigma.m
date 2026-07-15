function [sigma_log, sigma_alt_log] = VoI2sigma(VoI, D, t)
% Inverse: V/I -> log10 apparent conductivity for coincident loop
% Uses Newton-Raphson on the exact gamma function expression
%   V/I = (2/sigma) * gammainc(theta^2, 3/2)   where theta^2 = sigma * mu0 * D^2 / (4t)
%
% note that V(sigma) is NOT monotonic -- it peaks at theta^2 ~= 0.93. A given 
% V/I < V_peak corresponds to TWO half-space conductivities (LOW 
% and HIGH branch).
% Defaults to the LOW branch (standard TEM convention).
%
% D: loop side length (m)
% returns:
%   sigma_log:     log10 apparent conductivity (LOW branch)
%   sigma_alt_log: log10 conductivity on the OTHER branch (if different)

mu0 = 4 * pi * 10^(-7);

if VoI <= 0
    sigma_log = -6;
    if nargout > 1, sigma_alt_log = -6; end
    return;
end

% Late-time asymptotic (theta^2 << 1): V ~ sqrt(sigma)
coeff = 2 * (mu0 * D^2 / (4 * t))^(3/2) / 1.3293;
sigma_late = (VoI / coeff)^2;

% Early-time asymptotic (theta^2 >> 1): V ~= 2/sigma
sigma_early = 2 / VoI;

% Try Newton from BOTH initial guesses (handles non-monotonic V(sigma))
sigma_lo = sigma_late;
for k = 1:30
    t2 = sigma_lo * mu0 * D^2 / (4 * t);
    Vp = (2 / sigma_lo) * gammainc(t2, 3/2);
    d = Vp - VoI;
    if abs(d / max(VoI, 1e-30)) < 1e-10; break; end
    eps_s = max(sigma_lo * 1e-8, 1e-12);
    s2 = sigma_lo + eps_s;
    t2_2 = s2 * mu0 * D^2 / (4 * t);
    Vp2 = (2 / s2) * gammainc(t2_2, 3/2);
    dV = (Vp2 - Vp) / eps_s;
    if abs(dV) < 1e-30; break; end
    sigma_lo = sigma_lo - d / dV;
    sigma_lo = max(sigma_lo, 1e-12);
    sigma_lo = min(sigma_lo, 1e8);
end

sigma_hi = sigma_early;
for k = 1:30
    t2 = sigma_hi * mu0 * D^2 / (4 * t);
    Vp = (2 / sigma_hi) * gammainc(t2, 3/2);
    d = Vp - VoI;
    if abs(d / max(VoI, 1e-30)) < 1e-10; break; end
    eps_s = max(sigma_hi * 1e-8, 1e-12);
    s2 = sigma_hi + eps_s;
    t2_2 = s2 * mu0 * D^2 / (4 * t);
    Vp2 = (2 / s2) * gammainc(t2_2, 3/2);
    dV = (Vp2 - Vp) / eps_s;
    if abs(dV) < 1e-30; break; end
    sigma_hi = sigma_hi - d / dV;
    sigma_hi = max(sigma_hi, 1e-12);
    sigma_hi = min(sigma_hi, 1e8);
end

% If both converge to the same value (within 1% or near the peak), unique
if abs(sigma_lo - sigma_hi) / max(sigma_lo, sigma_hi) < 0.01
    sigma_log = log10(max(sigma_lo, 1e-12));
    if nargout > 1, sigma_alt_log = sigma_log; end
    return;
end

% Distinct branches: default to LOW branch (standard convention)
sigma_log = log10(max(sigma_lo, 1e-12));
if nargout > 1
    sigma_alt_log = log10(max(sigma_hi, 1e-12));
end
return