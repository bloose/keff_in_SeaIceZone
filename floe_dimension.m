function [L,d,P] = floe_dimension(sic,varargin),

% function [L,d,P] = floe_dimension(sic,varargin),
%
% floe_dimension: Calculate the characteristic floe length based on the
%                 the floe number distribution as function of sea ice 
%                  concentration.
%
% INPUT:  sic = Sea ice cover in %
%
% OUTPUT:
%   L = floe dimension (km)
%   d = logspace diameters
%   P = average floe area distribution (km2)
% 
% AUTHOR:  Brice Loose (osoberlice@gmail.com)
%
% REFERENCE:
%       Toyota, T., et al., (2006), "Characteristics of sea ice floe size 
%       distribution in the seasonal ice zone, Geophysical Research Letters"
%       33(L02616)L02616, doi: 10.1029/2005GL024556.
%
% DISCLAIMER:
%    This software is provided "as is" without warranty of any kind.  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%- Use the interpolate function to calculation the 
%- approximate floe dimension from the input Sea Ice Concentration (%).
%- from Toyota, we think the number distribution follows
%- N(d) per km2 = diameter^-alpha
%- alpha = 1.15 to 1.87
%- The cutoff for alpha = 



alpha1 = 1.7;
alpha2 = 1.15;

matchC = .04^(-alpha2+alpha1);
%d = (0.0001:1:10e3)/1e3;
d= logspace(log10(.005),log10(1),max(100,numel(sic)));
gt40 = d>.04;
lt40 = d<=.04;

%N = d.^-alpha1.*gt40;
%N2 = matchC.*d.^-alpha2.*lt40;
N =  matchC.*d.^-alpha1.*gt40;
N2 = d.^-alpha2.*lt40;


Nt = N+N2;

%- estimate circular and square surface area
A2 = (d/2).^2*pi;
A3 = d.*d;

%- These are like the Cumulative Distrib. Function between
%- diameter(d) as the abscissa and Flow Area Distribution as 
%- ordinate.
%SIC2 = cumsum(Nt.*A2.*gradient(d));
%SIC3 = cumsum(Nt.*A3.*gradient(d));

SIC2 = cumsum(Nt.*A2);
SIC3 = cumsum(Nt.*A3);


%mfloe = mean([SIC2; SIC3])*100;

%- These are like the PDF between
%- diameter(d) as the abscissa and Flow Area Distribution as 
%- ordinate.
%P = mean([Nt.*A2; Nt.*A3]);

P = Nt.*A2;

mfloe = cumsum(Nt.*A2)./sum(Nt.*A2);

L = interp1(mfloe*100,d,sic);

%- make a matrix where each column has the pmf for all values of x<=L

%P = interp1(d*1000,p,L);

fig = 'n';

if fig == 'y',
    figure(33); clf;
    plot(d*1000,matchC*(d).^-alpha1,'r');
    hold on;
    plot(d*1000,d.^-alpha2,'g');
    set(gca,'ylim',[.01 1e5]);
    set(gca,'xscale','log','yscale','log');
    set(gca,'xlim',[.1 1e4]);
    grid on;
    a = cumsum((d./2).^2*pi.*Nt)*diff(d(1:2));
    b = cumsum(d.^2.*Nt)*diff(d(1:2));
    figure;
    plot(mean([a;b]));
    mean([a;b])
end

return

