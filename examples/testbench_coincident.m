% simple testbench script, for 1D TEM (transient EM) occam inversion
% for coincident loop setup 
% DONG Hao
% 2010/01/07
% Yuxian, Hebei
% ======================================================================= %
clear
addpath(genpath('..'),'-end');

% some settings here
% terminating RMS misfit
Trms=1.5;
% number of maximum iteration
Niter=20; 
% loop diameter (loop side length, matching the original central-loop data)
D = 100;
% current 
I = 6;
% read a 19-layered model file
fid = fopen('coincide_loop.mod');
tmp = textscan(fid,'%f %f','CommentStyle', '#');
fclose(fid);
sigma0 = tmp{2};
z = tmp{1}(1:end-1);
% depth of each layer INTERFACE, note that z1=0
z0 = [0; z];
nz = length(z);
% read a TEM sounding data file containing apparent conductivity 
fid = fopen('coincide_loop.dat');
tmp = textscan(fid,'%f %f %f %f','CommentStyle', '#');
fclose(fid);
time = tmp{1};
asigma = tmp{2};
esigma = tmp{3}; 
VoI_obs = tmp{4};

% do the inversion
fprintf('=== Coincident loop Occam inversion ===\n');
[sigmai, asigmai, VoI_res] = occam1dtem(sigma0, z0, time, ...
    asigma, esigma, D, I, Trms, Niter, 'coincident');
% now try to read the true model 
fid = fopen('true.mod');
tmp = textscan(fid,'%f %f','CommentStyle', '#');
fclose(fid);
sigma_true = tmp{2};
z_true = tmp{1}(1:end-1);
z_true = [0; z_true];
% plot
figure(1);
plot1derr(time, 10.^asigma, VoI_obs, 'linear', 'rs', esigma);
hold on;
plot1derr(time, 10.^asigmai, VoI_res, 'linear', 'b-');
xlabel('Time (s)'); ylabel('V/I (V/A)');
title('TEM Response');
legend({'OBS', 'RESP'})
figure(2);
layers = diff(z0);
layers_true = diff(z_true);
plotlayer_log(sigma_true, layers_true, 'r-');
hold on;
plotlayer_log(sigmai, layers, 'b-');
xlabel('Conductivity (S/m)'); ylabel('Depth (m)');
title('1D Layered Model');
