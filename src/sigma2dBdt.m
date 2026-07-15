function dBdt = sigma2dBdt(sigma_a, D, I, t)
% Forward: log10 apparent conductivity -> dB/dt for central loop
% Uses the half-space central-loop formula
a = D / 2;
mu0 = 4 * pi * 10^(-7);
theta = sqrt((mu0 * sigma_a) ./ (4 * t));
x = theta * a;
dBdt = (I ./ (sigma_a * a^3)) .* (3 * erf(x) - (2/sqrt(pi)) .* x .* (3 + 2 * x.^2) .* exp(-x.^2));
return