function VoI = sigma2VoI(sigma_a, D, t)
% EXACT half-space coincident-loop V/I via incomplete gamma function
% Raiche & Spies (1981) equation (24):
%   V/I = (mu0 D^2 / 2t) * gammainc(theta^2, 3/2) / theta^2
%       = (2/sigma) * gammainc(theta^2, 3/2)
% where theta^2 = sigma * mu0 * D^2 / (4t)
%
% D: loop side length (m)
% This is the analytic integral expression, valid for ALL time ranges.

mu0 = 4 * pi * 10^(-7);
nt = length(t);
VoI = zeros(nt, 1);

for i = 1:nt
    sigma_i = sigma_a(i);
    t_i = t(i);
    
    if sigma_i <= 0 || t_i <= 0
        VoI(i) = 0;
        continue;
    end
    
    theta2 = sigma_i * mu0 * D^2 / (4 * t_i);
    VoI(i) = (2 / sigma_i) * gammainc(theta2, 3/2);
end
return