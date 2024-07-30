% simple testbench script, for 1D TEM (transient EM) occam inversion
% DONG Hao
% 2010/01/07
% Yuxian, Hebei
% ======================================================================= %
clear
addpath(genpath('..'),'-end');
% some settings here
% terminating RMS misfit
Trms=1.2;
% number of maximum iteration
Niter=20; 
% equivalent loop area
A = 100;
% current 
I = 6;
% read a 19-layered model file
fid = fopen('central_loop.mod');
tmp = textscan(fid,'%f %f','CommentStyle', '#');
fclose(fid);
sigma0 = tmp{2};
z = tmp{1}(1:end-1);
% depth of each layer INTERFACE, note that z1=0
z = [0; z];
nz = length(z);
% read a TEM sounding data file containing apparent conductivity 
fid = fopen('central_loop.dat');
tmp = textscan(fid,'%f %f %f %f','CommentStyle', '#');
fclose(fid);
time = tmp{1};
asigma = tmp{2};
esigma = tmp{3}; 
dBdt_obs = tmp{4};
% do the inversion
[sigmai,asigmai,dBdt_res]=occam1dtem(sigma0, z, time, asigma, esigma, A, I, Trms, Niter);
layers = diff(z);
% now try to read the true model 
fid = fopen('true.mod');
tmp = textscan(fid,'%f %f','CommentStyle', '#');
fclose(fid);
sigma_true = tmp{2};
z_true = tmp{1}(1:end-1);
z_true = [0; z_true];
layers_true = diff(z_true);
% ======================================================================= %
% now plot the results
figure(1);
plot1derr(time, 10.^asigma, dBdt_obs, 'linear', 'rs', esigma);
hold on 
plot1derr(time, 10.^asigmai, dBdt_res, 'linear', 'b-');
legend({'OBS', 'RESP'})
figure(2);
plotlayer_log(sigma_true, layers_true,'r-');
hold on
plotlayer_log(sigmai,layers,'b-');
% hasta la vista(?