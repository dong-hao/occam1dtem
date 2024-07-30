function [sigma,sigma_a,dBdt]=occam1dtem(sigma0, z, t, sigma_a0, sigma_e0,...
    A, I, Trms, Niter)
% main function of occam inversion
% for 1D TEM in a central loop set-up, note this inverts the apparent
% conductivity, instead of the H or dBdt
% in log10 space, i.e. both sigma and sigma_a are in log10 space
% DONG Hao
% 2010/01/07
% Yuxian, Hebei
% Note: this is but a toy scheme (no one is using 1D anymore, oh wait...) 
% when I was fiddling my 3D inversion code for my PhD. Thesis. 
% also note that I didn't use any sparse matrix here (need that if you are
% working on higher dimensions).
% 
% this is modified from my Occam MT routine
% Although this sort of worked (as far as I recall), I have not extensively 
% tested this script - please use with caution. 
%=========================================================================%
% input and output parametsigma:
% 
% sigma0:   array of conductivity of each layer of the starting model
%           in log10 space
% z:        array of layer DEPTH of each layer INTERFACE of the starting
%           model, z1=0(surface of the earth)
% t:        array of input data time step 
% a_sigma0: array of input data apparent conductivity
%           in log10 space
% e_sigma0: array of input data apparent conductivity (absolute) error 
% A:        (equivalent) loop area of the magnetic dipole source
% I:        electrical current of the source
% Trms:     target Root Mean Square misfit for iteration
% Niter:    maximum number of iteration
%=========================================================================%
% other parametsigma that might be useful:
%
% Nz:       number of model layers
% Nt:       number of frequencies 
% lambda0:  initial lagrange multiplier
%=========================================================================%
if nargin < 8
    Trms = 1.5;
    Niter = 20;
elseif nargin < 9
    Niter = 20;
end
% setup some params
M=length(sigma0);
lambda=10;
dlambda=10;
% initial sigma
sigma = sigma0;
%==================initialize inversion iteration=========================%
disp('========== starting Occam iterations ==========')
for iter=1:Niter
    % the outer loop
    disp(['================ iteration ', num2str(iter),' =================='])
%   plot_along(sigma,z,iter);
    lambdar=lambda*dlambda; % lambda on the right
    lambdal=lambda/dlambda; % lambda on the left
    % first get the Jacobian and the 
    [J,sigma_a]=jacob10(sigma,z,t,A,I);
    imrms=rms1(sigma_a,sigma_a0,sigma_e0);     % starting rms
    fprintf('! previous RMS = %5.3f \n', imrms);
    dsigma=sigma_a0-sigma_a;
    % one-step occam iteration for three different lambdas 
    sigmam=occam(dsigma,sigma_e0,sigma,lambda, J); % conductivity in the middle
    sigmar=occam(dsigma,sigma_e0,sigma,lambdar,J); % conductivity on the right
    sigmal=occam(dsigma,sigma_e0,sigma,lambdal,J); % conductivity on the left
    sigma_am = tem1dfwd10(sigmam, z, t, A, I);
    sigma_ar = tem1dfwd10(sigmar, z, t, A, I);
    sigma_al = tem1dfwd10(sigmal, z, t, A, I);
    [ChiSl,ChiSm,ChiSr]=dispfit(sigma_a0, sigma_am, sigma_al, ...
        sigma_ar,sigma_e0, sigmal ,sigmam, sigmar, lambdal ,lambda, lambdar);  
    for ifind=1:10
        % the inner loop to sweep through different Lagrange multipliers
        % (or mu in origin Occam theory)
        if ChiSm <= ChiSl && ChiSm <= ChiSr
            disp('============regional minimum found=============')
            break
        elseif ChiSm > ChiSl && ChiSr > ChiSl
            disp('================<< searching <<================')
            if lambdal/dlambda<=0.0001
                disp('! lambda is too small, restart iteration')
                break
            end 
            sigmar=sigmam;
            lambdar=lambda;
            sigmam=sigmal;
            lambda=lambdal;
            lambdal=lambda/dlambda;
            sigmal=occam(dsigma,sigma_e0,sigma,lambdal,J);  % conductivity on the left
        elseif ChiSr < ChiSm && ChiSr < ChiSl
            disp('================>> searching >>================')
            if lambda*dlambda>=100000
                disp('! lambda is too large, restart iteration')
                break
            end
            sigmal=sigmam;
            lambdal=lambda;
            sigmam=sigmar;
            lambda=lambdar;
            lambdar=lambda*dlambda;
            sigmar=occam(dsigma,sigma_e0,sigma,lambdar,J);  % sigmaistivity on the right
        else 
            % seems not converging
            % use the middle value for next inversion iteration
            disp('=========cannot find a local minimum===========')
            break
        end
        sigma_am = tem1dfwd10(sigmam,z,t,A,I);
        sigma_al = tem1dfwd10(sigmal,z,t,A,I);
        sigma_ar = tem1dfwd10(sigmar,z,t,A,I);
        [ChiSl,ChiSm,ChiSr]=dispfit(sigma_a0, sigma_am, sigma_al, ...
        sigma_ar,sigma_e0, sigmal ,sigmam, sigmar, lambdal ,lambda, lambdar);   
    end
    sigma=sigmam;
    sigma_a=sigma_am;
    irms=rms1(sigma_a,sigma_a0,sigma_e0);
    if irms<=imrms
        % So far so good...
        disp(['! finishing iteration # ' num2str(iter)]);
    else
        % Huston, we have a convergence problem. 
        disp(['! warning: iteration # ' num2str(iter) ' not converged']);
        disp('! try stablizing inversion with smaller lambda...');
        disp(['! RMS= ' num2str(irms)]);
        lambda=lambda/dlambda.^2;
        continue
    end
    if irms<=Trms
        % check if the desired rms is reached
        fprintf('! current RMS = %5.3f \n', irms);
        fprintf('! target rms (%5.3f given by user) reached \n', Trms);
        disp('! try to find a smoothest model that fit the data')
        rmsr=irms;
        dsigma=sigma_a0-sigma_a;
        while rmsr<=Trms && lambda < 100000
            disp('================>> searching >>================')
            lambda=lambda*sqrt(dlambda);
            sigmar=occam(dsigma,sigma_e0,sigma,lambda,J);
            sigma_ar = tem1dfwd10(sigmar,z,t,A,I);
            rmsr=rms1(sigma_a0,sigma_ar,sigma_e0);
            fprintf('! evaluating a smoother model with Lambda = %5.3f \n', lambda);
            fprintf('! current RMS = %5.3f \n', rmsr);
            if rmsr < Trms 
                fprintf('! model accepted \n');
                sigma_am = sigma_ar;
                sigmam = sigmar;
            else
                fprintf('! model rejected \n');
            end
        end
        disp('! picking up the smoothest model...')
        disp('! exiting...')
        sigma=sigmam;        
        break
    end
    if ifind>=10
        disp('! maximum lambda search reached, start next search...')
    end
end
% =========== end main iterations here =========== %
if iter>=Niter
    disp('! iteration limit reached, stop...')
else
    sigma_a = sigma_am;
end
if nargout > 2 
    % exciplitly output a dBdt value
    [sigma_a, dBdt] = tem1dfwd10(sigma, z, t, A, I);
end
frms=rms1(sigma_a,sigma_a0,sigma_e0);  
ChiS=chi2(sigma_a,sigma_a0,sigma_e0);
if ChiS<=2*M
    disp('! WARNING: desired chi2 reached. ');
    disp('! You might be over fitting the data! ')
end
fprintf('! final RMS = %5.3f \n', frms);
return

function sigma=occam(dpara,epara,sigma,lambda,J)
% main function for a single step of occam iteration 
% the dpara and epara be the diff and variance of the parameter
% here in 
N=size(J,2);% number of model
D=mk_pmat(N); % roughness matrix
% penalty matrix
W=diag(1./epara);
a=lambda*(D)'*D;%+diag(ones(N,1));
% construct weighted Jacobian premultiplied by TRANS(W.J)
wjtwj=(W*J)'*(W*J);
% construct weighted translated data premultiplied by TRANS(W.J)
wjtwd=(W*J)'*(W*(dpara+J*sigma));
sigma=(a+wjtwj)\wjtwd;
return

function DEL=mk_pmat(N)
% constructing roughness matrix 
% note this is for L2 norm
DEL=diag(ones(N,1),0)*2-diag(ones(N-1,1),-1)-diag(ones(N-1,1),1); 
DEL(1,1)=1;
DEL(N,N)=1;
return 

function [ChiSl,ChiSm,ChiSr]=dispfit(co, cm, cl, cr, stderr, sigmal, sigmam, sigmar, lambdal ,lambda, lambdar)
ChiSm=chi2(co,cm,stderr);
ChiSl=chi2(co,cl,stderr);
ChiSr=chi2(co,cr,stderr);
roughl=roughness1(sigmal,1);
roughm=roughness1(sigmam,1);
roughr=roughness1(sigmar,1);
fprintf('               L          M          R\n');
fprintf('Chi Square   : %8.4e %8.4e %8.4e \n',...
    ChiSl, ChiSm, ChiSr);
fprintf('Roughness(L1): %8.4e %8.4e %8.4e \n',...
    roughl, roughm, roughr);
fprintf('Lambda       : %8.5g %8.5g %8.5g \n',...
    lambdal, lambda, lambdar );
return