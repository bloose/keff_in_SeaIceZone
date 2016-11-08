%- Example script for providing inputs to keff_SIZ.m for calculating
%- the gas transfer velocity in m/d from shear and buoyancy driven
%- turbulence in the Ice-Ocean Boundary Layer.
%
% AUTHOR:  Brice Loose (osoberlice@gmail.com)
%          Arash Bigdeli (arash.big@gmail.com)
%
% REFERENCE:
%       Loose et al., (2014), "A parameter model of gas exchange for the 
%       seasonal sea ice zone", Ocean Ocean Sci., 10, 17-28, 2014
%       http://www.ocean-sci.net/10/17/2014/
%       doi:10.5194/os-10-17-2014
%
% DISCLAIMER:
%    This software is provided "as is" without warranty of any kind.  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% example 1 
clc;clear,close all;fig1 = figure;

%This example shows the relation of sea ice cover and speed vs gas exchange
Uwind = [5 10 15];      %Wind speed in m/s
Uice = 0.02.*Uwind;     % free drift 0.02*Uwind
waterT = 0;             %Water Temp in deg C - 1.5
airT   = 5;             %Air Temp in deg C 1
sic = 1:1:100;          %Sea ice Concentration in %.

for i = 1:length(Uwind)
[keff_i(:,i)] = keff_SIZ(Uice(i),Uwind(i),sic,waterT,airT);
[keff(:,i)] = keff_SIZ(0,Uwind(i),sic,waterT,airT);
end

%Plotting
ax1 = plot(sic,keff_i,'linewidth',2);hold on;
ax2 = plot(sic,keff_i,'linewidth',2);
ax3=plot(sic,keff,'--','linewidth',2);
ylabel('k_{eff} (m/d)');  xlabel('SI cover %');
title('Gas exchange vs SI cover')
ah=axes('position',get(gca,'position'),'visible','off');
legend(ax1,'Uwind 5ms^{-1}','Uwind 10ms^{-1}','Uwind 15ms^{-1}'...
    ,2,'location','northeast') 
legend(ah,[ax2(1) ax3(1)],'SI in free drift','Stationary SI'...
    ,2,'location','southwest')


