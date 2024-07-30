function plot1derr(time,var1,var2,opt,ccode,err1,err2)
% a barn door function to plot 1d apparent conductivity and dBdt...
% 
if nargin < 3
        error('not enough input arguments, 4 at least')  
elseif nargin < 4
	opt='log10';
	ccode='b-';
elseif nargin < 5
        ccode='b-';
end
if strcmp(opt,'linear')==1
    sig = log10(var1);
    dBdt = var2;  
    if nargin > 5
        sige = 0.434*err1;
    elseif nargin > 6
        sige = 0.434*err1;
        dBdte = 0.434*err2;
    end 
else % in log10 domain
    sig = var1;
    dBdt = var2;
    if nargin > 5
        sige = err1;
    elseif nargin > 6
        sige = err1;
        dBdte = err2;
    end
end
a1=subplot(2,1,1);
if nargin <= 5
    plot(a1,time,sig,ccode);
else
    errorbar(a1,time,sig,sige,ccode);
end
hold(a1,'on');
set(a1,'xscale','log');
set(a1,'xgrid','on','ygrid','on')
xlabel('time(s)');
ylabel('log10 app. conductivity (S/m)')
a2=subplot(2,1,2);
if nargin <= 6
    plot(a2,time,dBdt,ccode);
else
    errorbar(a2,time,dBdt,dBdte,ccode);
end
hold(a2,'on');
set(a2,'xscale','log');
set(a2,'xgrid','on','ygrid','on')
set(a2,'yscale','log');
xlabel('time(s)')
ylabel('dB/dt (A/m^2)')
return