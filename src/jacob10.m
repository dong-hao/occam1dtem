function [J, sigma_a] = jacob10(sigma, z, time, D, I)
% Analytical Jacobian for 1D TEM central loop (ABFM forward model)
% used for occam 1D TEM inversion
% Uses implicit differentiation of the ABFM fixed point
% instead of the finite difference approximation in the previous 
% version.
% 
% specifically, we have:
% J(i,j) = d log10(sigma_a_i) / d log10(sigma_j)
%          (i = time channel, j = layer)
%
% DONG Hao
% 2010/1/7
% Yuxian, Hebei
%
NL = length(sigma);
NT = length(time);
J = zeros(NT, NL);

% some duplicated code from tem1dfwd10.m
mu0 = 4 * pi * 1e-7;
alpha = 0.6;
c_const = 1.2;
z_pad = [z; Inf];

sigma_lin = 10.^sigma;

% ABFM forward with derivative - can be simplified as we already have the 
% forward calculation in tem1dfwd10.m 
sigma_a_lin = mean(sigma_lin) * ones(NT, 1);
relres = 1;

k = 1;
while (relres > 1e-6 && k < 30)
    d = sqrt((c_const * time) ./ (sigma_a_lin * mu0));
    F = zeros(NT, NL);
    dF_dd = zeros(NT, NL);

    for i = 1:NT % loop over time channels
        for j = 1:NL % loop over layers
            z_top = z_pad(j);
            z_bot = z_pad(j+1);
            % F1/dF1 at the top interface
            if z_top <= d(i)
                uz = z_top / d(i);
                F1 = (2 - uz) * uz;
                dF1_dd = -2*z_top/d(i)^2 + 2*z_top^2/d(i)^3;
            else
                F1 = 1;
                dF1_dd = 0;
            end
            % F2/dF2 at the bottom interface
            if isinf(z_bot)
                F2 = 1;
                dF2_dd = 0;
            elseif z_bot <= d(i)
                uz = z_bot / d(i);
                F2 = (2 - uz) * uz;
                dF2_dd = -2*z_bot/d(i)^2 + 2*z_bot^2/d(i)^3;
            else
                F2 = 1;
                dF2_dd = 0;
            end

            F(i, j) = F2 - F1;
            dF_dd(i, j) = dF2_dd - dF1_dd;
        end
    end

    app_sigma1 = sigma_lin .* F';
    av_sigma_old = sigma_a_lin;
    sigma_a_lin = sum(app_sigma1, 1)';
    sigma_a_lin = alpha * sigma_a_lin + (1 - alpha) * av_sigma_old;
    relres = norm(av_sigma_old - sigma_a_lin) / norm(av_sigma_old);
    k = k + 1;
end

% ---------- Analytical Jacobian via implicit differentiation ---------- %
% At the fixed point we have sigma_a = sum_j ( sigma_j * F_j(sigma_a) )
% Implicit differentiation gives:
%   d sigma_a / d sigma_j = F_j / ( 1 - sum_k sigma_k * dF_k/d(sigma_a) )
%
% Chain rule: 
%   dF_k/d(sigma_a) = dF_k/dd * dd/d(sigma_a)
%   d = sqrt(c*t / (mu0 * sigma_a))  -->  dd/d(sigma_a) = -d / (2 * sigma_a)
%
% So we have:  
%       C = 1 - sum_k ( sigma_k * dF_k/dd * (-d/(2*sigma_a)) )
%         = 1 + (d/(2*sigma_a)) * sum_k ( sigma_k * dF_k/dd )

for i = 1:NT
    d_i = d(i);
    s_a_i = sigma_a_lin(i);
    C_i = 1 + (d_i / (2 * s_a_i)) * sum(sigma_lin .* dF_dd(i, :)');
    % J(i,j) = (sigma_j / sigma_a_i) * (F_ij / C_i)
    J(i, :) = (sigma_lin' / s_a_i) .* (F(i, :) / C_i);
end

sigma_a = log10(sigma_a_lin);
end
