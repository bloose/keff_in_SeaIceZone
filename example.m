%- Example script for providing inputs to keff_SIZ.m for calculating
%- the gas transfer velocity in m/d from shear and buoyancy driven
%- turbulence in the Ice-Ocean Boundary Layer.
%
% AUTHOR:  Brice Loose (osoberlice@gmail.com)
%
% REFERENCE:
%       Loose et al., (2013), "A parameter model of gas exchange for the 
%       seasonal sea ice zone", Ocean Sci. Discuss., 10(4), 1169-1204.
%
% DISCLAIMER:
%    This software is provided "as is" without warranty of any kind.  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Uwind = 8;     %Wind speed in m/s
Uice = .08;    %Ice velocity in m/s
waterT = -1.5;  %Water Temp in deg C
airT   = 1;  %Air Temp in deg C
sic = 1:1:100;  %Sea ice Concentration in %.

[keff_ice, keff_wave] = keff_SIZ(Uice,Uwind,sic,waterT,airT);