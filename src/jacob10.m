function [J, sigma_a] = jacob10(sigma, z, time, A, I)
% calculate the Jacobian dsigma_a/dsigma for the 1D TEM case
% with the hard way (with an arbitary pertubation)
% for dsigma_a/dsigma
% used for occam 1D TEM inversion
% DONG Hao
% 2010/1/7
% Yuxian, Hebei
% fixed permeability, for now
% note that the Jacobian is only related to conductivity (i.e. useful for 
% smooth inverions only)
NL = length(sigma);
NT = length(time);
% a small value for pertubation
pert = 1e-6;
J = zeros(NT, NL);
sigma_a0 = tem1dfwd10(sigma, z, time, A, I);
for i = 1:NL % m forward modellings 
    psigma = sigma;
    psigma(i) = (1 + pert) * sigma(i);
    sigma_a = tem1dfwd10(psigma, z, time, A, I);
    psigma_a = sigma_a - sigma_a0;
    % dsigma_a/dsigma
    J(:, i) = psigma_a / (pert * sigma(i));
end