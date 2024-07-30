function plotlayer_log(res,layer,ccode)
% a barn door function to plot 1d layered model for 1D TEM inversion
% in logrithm scale
if nargin<3
    ccode='k';
end
if (size(res,2)~=size(layer,2))
    if length(res)==length(layer)
        layer=layer';
    else
	    disp('please check the input parametres ')
        error('res and layer should have the same size');
    end
end
NL=length(layer);
depth=0.1;
h=gca;
for i=1:NL
    plot(h,[res(i) res(i)],[depth depth+layer(i)],ccode,'linewidth',2)
    hold(h,'on');
    plot(h,[res(i) res(i+1)],[depth+layer(i) depth+layer(i)],ccode,'linewidth',2)
    depth=depth+layer(i);
end
plot(h,[res(NL+1) res(NL+1)],[depth depth+layer(NL)*10],ccode,'linewidth',2);
hold(h,'off');
set(h, 'ydir', 'reverse');
set(h, 'xlim', [-4 0])
set(h, 'ylim', [1 depth+layer(NL)*10]);
xlabel(h,'log_{10} Conductivity(S/m)');
ylabel(h,'Depth(m)');
set(h,'yscale','log');
grid on;
