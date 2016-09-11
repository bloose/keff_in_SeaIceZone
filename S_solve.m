function b=S_solve(a)

% function [b] = S_solve(a,varargin),
%
% S_solve: Calculations for the heat and buoyancy flux to/from the ice.
% Based on interface fluxes using the 3 equation solution (IOBL, 6.9). 
% <w'T'>_0 = w_0*Q_L + q
% <w'S'>_0 = w_0*(Sice-S0)
%
% INPUT:
% a.us0 = Friction velocity (m/s)
% a.alpha_h = Thermal expancion coefficient (T^-1)
% a.alpha_S = Haline contraction coefficient (S^-1)
% a.Tml = Mixed-layer temperature (deg C);
% a.Sml = Mixed-layer salinity (psu);
% a.Sice = Ice salinity
% a.wperc = Percolation velocity (m/s);
%
% OUTPUT: b same as a except add b.T0,b.S0,b.w0,b.QL,b.m.b.rho_cp
%         b.wb0 = Buoyancy flux in/out of ice in m^2/s^3 or W/kg
% 
% AUTHOR:  Miles McPhee
%
% REFERENCE:
%       McPhee (2010) "Air-ice-Ocean Interaction: Turbulent Ocean
%       Boundary Layer Exchange Processes", Springer.
%
% DISCLAIMER:
%    This software is provided "as is" without warranty of any kind.  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% b=S_solve(a);
% 26 Feb 09 MGM 
% INPUT a.us0
%         .Tml
%         .Sml
%         .Sice
%         .Hup_ice    upward heat flux in the ice, W m^-2
%         .alpha_h
%         .alpha_S
%         .wperc
% OUTPUT: b same as a except add b.T0,b.S0,b.w0,b.QL,b.m.b.rho_cp
Lfresh=333.5e3;
gee=9.8;
%cp=specific_heat(a.Tml,a.Sml,0);
cp=gsw_cp_t_exact(a.Sml,a.Tml,0);
rho=gsw_rho_t_exact(a.Tml,a.Sml,0);
a.rho_cp=rho.*cp;
a.QL=(1-0.03*a.Sice).*Lfresh./cp;
qdot=a.Hup_ice./a.rho_cp;
a.m=-gsw_t_freezing(a.Sml,0)./a.Sml;
T_H=a.Tml-qdot./a.alpha_h./a.us0;
T_L=a.alpha_S.*a.QL./a.alpha_h;
T_p=a.wperc.*a.QL./a.alpha_h./a.us0;
A=a.m;
B=T_H+T_L+T_p-a.m.*a.Sice;
C=-(T_H+T_p).*a.Sice-T_L.*a.Sml;
S0=(-B+sqrt(B.^2-4.*A.*C))./2./A;
w0=(a.alpha_h.*a.us0.*(a.Tml+a.m.*S0)-qdot)./a.QL;
b=a;
b.S0=S0;
b.T0=-b.m.*S0;
b.w0=w0;
b.wt0=b.w0.*b.QL+b.Hup_ice./b.rho_cp;
b.ws0=b.w0.*(b.S0-b.Sice);
%[bt,bs]=beta_ts(b.T0,b.S0);
bs = gsw_beta_const_t_exact(b.S0,b.T0,0);
bt = gsw_alpha_wrt_t_exact(b.S0,b.T0,0);
b.wb0=gee.*(-bt.*b.wt0+bs.*b.ws0);