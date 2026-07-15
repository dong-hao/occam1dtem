function sigma_log = dBdt2sigma(dBdt, D, I, t)
% Inverse: dB/dt -> log10 apparent conductivity for central loop
% Uses Newton-Raphson on the half-space central-loop formula
%   dB/dt = (I/(sigma*a^3)) * (3*erf(x) - (2/sqrt(pi))*x*(3+2*x^2)*exp(-x^2))
%   where x = a*sqrt(mu0*sigma/(4*t)), a = D/2
%
% The central-loop function is monotonic in sigma (unlike coincident loop),
% so the inverse is unique for each time channel.
%
% D: transmitter loop DIAMETER (m)
% I: transmitter current (A)
% t: time (s), may be a vector
% returns:
%   sigma_log: log10 apparent conductivity (same size as t)

mu0 = 4 * pi * 10^(-7);
a = D / 2;
nt = length(t);
sigma_log = zeros(nt, 1);

for i = 1:nt
    dB_i = dBdt(i);
    t_i = t(i);

    if dB_i <= 0
        sigma_log(i) = -6;
        continue;
    end

    % Late-time asymptotic (x << 1): g(x) ~ (8/5√π)*x^5
    %   dB/dt = I*a^2 * (8/(5√π)) * (μ₀/(4t))^(5/2) * σ^(3/2)
    K = I * a^2 * (8/(5*sqrt(pi))) * (mu0/(4*t_i))^(5/2);
    sigma_late = (dB_i / K)^(2/3);

    % Early-time asymptotic (x >> 1): g(x) → 3, dB/dt ≈ 3I/(σ*a^3)
    sigma_early = 3 * I / (a^3 * dB_i);

    % Newton-Raphson from the late-time initial guess
    sigma_lo = sigma_late;
    for k = 1:30
        x = a * sqrt(mu0 * sigma_lo / (4 * t_i));
        g = 3 * erf(x) - (2/sqrt(pi)) * x * (3 + 2 * x^2) * exp(-x^2);
        Vp = (I / (sigma_lo * a^3)) * g;
        d = Vp - dB_i;
        if abs(d / max(dB_i, 1e-30)) < 1e-10; break; end
        eps_s = max(sigma_lo * 1e-8, 1e-12);
        s2 = sigma_lo + eps_s;
        x2 = a * sqrt(mu0 * s2 / (4 * t_i));
        g2 = 3 * erf(x2) - (2/sqrt(pi)) * x2 * (3 + 2 * x2^2) * exp(-x2^2);
        Vp2 = (I / (s2 * a^3)) * g2;
        dV = (Vp2 - Vp) / eps_s;
        if abs(dV) < 1e-30; break; end
        sigma_lo = sigma_lo - d / dV;
        sigma_lo = max(sigma_lo, 1e-12);
        sigma_lo = min(sigma_lo, 1e8);
    end

    % Newton-Raphson from the early-time initial guess (backup)
    sigma_hi = sigma_early;
    for k = 1:30
        x = a * sqrt(mu0 * sigma_hi / (4 * t_i));
        g = 3 * erf(x) - (2/sqrt(pi)) * x * (3 + 2 * x^2) * exp(-x^2);
        Vp = (I / (sigma_hi * a^3)) * g;
        d = Vp - dB_i;
        if abs(d / max(dB_i, 1e-30)) < 1e-10; break; end
        eps_s = max(sigma_hi * 1e-8, 1e-12);
        s2 = sigma_hi + eps_s;
        x2 = a * sqrt(mu0 * s2 / (4 * t_i));
        g2 = 3 * erf(x2) - (2/sqrt(pi)) * x2 * (3 + 2 * x2^2) * exp(-x2^2);
        Vp2 = (I / (s2 * a^3)) * g2;
        dV = (Vp2 - Vp) / eps_s;
        if abs(dV) < 1e-30; break; end
        sigma_hi = sigma_hi - d / dV;
        sigma_hi = max(sigma_hi, 1e-12);
        sigma_hi = min(sigma_hi, 1e8);
    end

    % Both converge to the same value since central-loop is monotonic
    sigma_log(i) = log10(max((sigma_lo + sigma_hi) / 2, 1e-12));
end
return
