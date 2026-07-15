function [sigma_a, out2] = tem1dfwd10(sigma, z, t, D, I, config)
% a barn door 1d layered forward routine for 
% TEM 1D inversion
% in log10 space, i.e. both sigma and sigma_a are in log10 space
% 
% DONG Hao
% 2010/1/7
% Yuxian, Hebei
%
% sigma -> conductivity (nlayer by 1 array)
% z     -> layer interface depth z(0) = 0 for ground-based cases 
%           (nlayer by 1 array)
% t     -> time step
% D     -> central-loop mode: transmitter loop DIAMETER (m)
%          coincident-loop mode: loop side length (m)
% I     -> electrical current of the source
% config-> 0 or 'central'   (default): central loop, out2 = dB/dt
%          1 or 'coincident'        : coincident loop, out2 = V/I

% See: Christensen, N. B. A Generic 1-D Imaging Method for Transient
% Electromagnetic Data. Geophysics. 2002. 67, 438-447.
% this is essentially an approximate method which express a layered
% response with a series of half-space responses

% check data consistency
if nargin < 6
    config = 0; % default to central loop
end
if ischar(config)
    if strcmpi(config, 'coincident')
        config = 1;
    else
        config = 0;
    end
end

if size(sigma,1) ~= size(z,1)
    error('inconsistence size of resistivity and layer number...')
else
    z = [z; Inf];
end
nlayer = size(z,1);

% fixed permeability (air), for now
mu0 = 4 * pi * 10^(-7); 

% factors
alpha = 0.6;
c = 1.2;

% convert to (initial) true conductivity 
sigma_a = zeros(nlayer, 1);
sigma = 10.^sigma;

% initial average sigma for each time step
for i = 1 : length(t)
    sigma_a(i) = mean(sigma);
end

% some initial setup
relres = 1;
k = 1;
while (relres > 1e-6 && k < 30) % hard coded here...
    d = sqrt((c*t)./(sigma_a*mu0));
    F1 = zeros(length(t),(length(z)-1));
    F2 = zeros(length(t),(length(z)-1));
    F  = zeros(length(t),length(z)-1);
    for i = 1 : length(t)
        for j = 1 : (length(z)-1)
            if z(j) <= d(i) 
                F1(i,j) = (2-(z(j)/d(i))).*(z(j)/d(i));
            elseif z(j)> d(i)
                F1(i,j) = 1;
            end
           if z(j+1) <= d(i) 
                F2(i,j) = (2-(z(j+1)/d(i))).*(z(j+1)/d(i));
            elseif z(j+1)> d(i)
                F2(i,j) = 1;
            end
            F(i,j)=F2(i,j)-F1(i,j);    
        end
    end
    app_sigma1 = sigma.*F';
    av_sigma_old  = sigma_a;
    sigma_a   = sum(app_sigma1,1)';
    sigma_a   = alpha * sigma_a + (1 - alpha) * av_sigma_old;
    relres = norm(av_sigma_old - sigma_a)./norm(av_sigma_old);
    k = k + 1;
end

if nargout > 1
    if config == 0
        % central loop: compute dB/dt from half-space formula
        out2 = sigma2dBdt(sigma_a, D, I, t);
    else
        % coincident loop: compute V/I from half-space formula
        % D is the loop side length (Raiche & Spies, 1981)
        out2 = sigma2VoI(sigma_a, D, t);
    end
end
% convert back to log10 space
sigma_a = log10(sigma_a);
return

