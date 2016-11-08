function [keff,eps] = keff_SIZ(Uice,Uwind,sic,waterT,airT,varargin),

% function [keff,eps] = keff_SIZ_swh(Uice,Uwind,sic,waterT,airT,varargin),
%
% keff_SIZ: The effective gas transfer velocity in the sea ice zone.
%
% INPUT:
%       Uice = ice velocity (m/s)
%       Uwind = wind speed (not vector) in m/s
%       sic = sea ice concentration (%)
%       waterT = water temp in deg C
%       airT = air temp in deg C.
%
% OUTPUT:
%   keff = effective gas transfer velocity from sea ice-driven turbulence  [m/d]
%   s2keff = effective gas transfer velocity from capillary-gravity waves  [m/d]
%
% DEPENDENCIES:
%   These routines require the Gibbs Seawater toolbox:
%   http://www.teos-10.org/
%
% AUTHOR:  Brice Loose (osoberlice@gmail.com)
%
% REFERENCE:
%       Loose et al., (2014), "A parameter model of gas exchange for the
%       seasonal sea ice zone", Ocean Sci., 10(4), doi:10.5194/os-10-1-2014.
%
% DISCLAIMER:
%    This software is provided "as is" without warranty of any kind.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if nargin > 6,
    hum = varargin{1};
else
    hum = 273.16./8.3141./satvap(waterT(:))./1e2*.06;
end
if nargin > 6,
    %- Set the Mixed-layer depth (m) if it has not been specified.
    mld = varargin{2};
else
    mld = 30;
end
if nargin > 8
    %- Set the salinity if it has not been specified.
    waterS = varargin{3};
else
    waterS = repmat(32,size(waterT));
end

%- SIC (as decimal)
Ai = sic(:)./100;
airT = airT(:);
waterT = waterT(:);
waterS = waterS(:);
Uice = Uice(:);
Uwind = Uwind(:);

%- convert to Kelvin
waterT = waterT+273.16;
airT = airT+273.16;

%--- COEFICIENTS
g = 9.8;
cd = 1.5e-3;    %- Air-water drag coefficient
vis = 1.787e-6; %- Kinematic viscosity m2/s
Sc = 660;       %- Schmidt Number (unitless)
n = 0.5;        %- Schmidt Number exponent
% z = 3;
z = 3;
z0 = .01;       %- Roughness length
acp  = 1.005e3; %- Specific heat of air (J/kg/K)
Lh = 2.26e6;    %- Latent heat of melting (J/kg)
rho_a = 1.013e5/287.058./(airT);
rho_w = gsw_rho_t_exact(waterS, waterT-273.16,0);




%% --------------------------------------------------------------------- %%
%----------------- (1) BUOYANCY FLUX CALCULATION -------------------------%
%-------------------------------------------------------------------------%

%- J_b^o = g*Ch*U*(Tw-T)/T -- (m2/s3)


%- heat capacity and expansion/contraction coefs.
%wcp  = sw_cp(waterS,waterT-273.16,0);
wcp = gsw_cp_t_exact(waterS,waterT-273.16,0);
%alpha = sw_alpha(waterS,waterT-273.16,0);
alpha = gsw_alpha_wrt_t_exact(waterS,waterT-273.16,0);
%beta = sw_beta(waterS,waterT-273.16,0);
beta = gsw_beta_const_t_exact(waterS,waterT-273.16,0);

%- precalculate u*0 because we need for ice-water
us = Uice.*.41/log(z/z0);
sens_hf = rho_a*acp*cd.*abs(Uwind).*(waterT-airT);
Jb_sens = 9.8*alpha.*sens_hf./rho_a./acp;



%--- Andreas and Murphy (1986), "Bulk transfer coefficients for heat and
%- momentum over leads and polynyas". JPO, V. 16, 1875-1883.
%- Qs = Qsat(ts).  Then Qr is "ref" specific humidity at same height as
%- Uwind.
%- Units are kg-w/m3-air
%- pv = nrt --- n/v(mol/m3)=rt/p.  n/v*.06kg/mol
Qs = 273.16./8.3141./satvap(waterT-273.15)./1e2*.06;

%- Buoyancy flux due to bulk evaporative Latent heat flux.
JqLf = Lh.*cd.*abs(Uwind).*(Qs-hum).*rho_a;
Jb_Lf = g.*alpha.*JqLf./rho_a./acp;

%- Buoyancy flux due to salinification from evaporation
Jb_salt = g.*beta.*waterS.*(JqLf./Lh)./1e3;


%- Buoyancy flux from ice melt/freeze;
%- Based on interface fluxes using the 3 equation solution (IOBL, 6.9).
%- <w'T'>_0 = w_0*Q_L + q
%- <w'S'>_0 = w_0*(Sice-S0)
%-
%- Reference: McPhee (2010) Air-ice-Ocean Interaction, Springer.
a.us0 = us;
a.alpha_h = .0057;
a.alpha_S = 4e-4;
a.Tml = waterT-273.16;
a.Sml = waterS;
a.Sice = 6;     %- ice salinity
a.wperc = 1e-7;

%- calculate ice temperature as half way between water and air temp.
Tice = nanmean([gsw_t_freezing(waterS,0)+273.16 airT],2);
%Tice = nanmean([gsw_t_freezing(waterS,0)+273.16 airT],2);
%- Kice ~ Kfresh + beta*Sice/Tice
Kice = 2.04 + 0.117*a.Sice./Tice;

a.Hup_ice = -Kice.*(Tice-waterT);

b = S_solve(a);

%- Output b.wb0 is already in W/kg or m2/s3.
Jb_ice = b.wb0;


%- Positive Flux out of the ocean, corresponds to a loss of buoyancy
%- This is why first terms have "negative" signs in front.
Jb = nansum([-(1-Ai).*Jb_sens,-(1-Ai).*Jb_Lf, (1-Ai).*Jb_salt, +Ai.*Jb_ice],2);


%% --------------------------------------------------------------------- %%
%----------------- (2) SHEAR CALCULATION ---------------------------------%
%-------------------------------------------------------------------------%

%- CALC tau_skin-iw the ice/water skin friction using the Rossby similarity
%- theorem:  Vo/u* = 1/k(log(u*/fzo) - A + iB)
%--------------------------------------------
%- Reference: McPhee (2010) Air-ice-Ocean Interaction, Springer.

rho_a = 1.013e5/287.058./(airT);
tau_aw = (1-.14)*rho_a.*cd.*(Uwind.^2);
us_aw = sqrt(tau_aw./rho_w);
% now set roughness to 1 cm, Coriolis parameter for 75N
kappa = .41;
f=1.4e-4;
Rostar=b.us0./f./z0;

% compute the Obukhov length
%L=b.us0.^3./kappa./b.wb0;
L=b.us0.^3./kappa./Jb;
Lo = L;
L(L<1) = 1;

L = 1;

% mu* is ratio of planetary to Obukhov scales
mustar=b.us0./f./L;

%- When mustar, stability parameter is negative, stratification is not
%- affecting surface shear.
mustar(mustar<0)=0;

[Und,etastar,A,B]=U0(Rostar,mustar);

Vstable=b.us0.*Und./etastar;
Vsr = real(Vstable).*kappa./(log(Rostar)-A);
Vsi = imag(Vstable).*kappa./(log(Rostar)-A-B);
us_iw = sqrt(Vsr.^2+Vsi.^2);
tau_iw = us_iw.^2.*rho_w;


%- CALC FORM DRAG FROM ESTIMATE OF SEA ICE COVER AND DRIFT VELOCITY -%
%-- This takes into account the recognition that the flow size
%- distribution is a pmf w.r.t. sea ice concentration.

%- Calculate L from SIC:
[floe_dim, dia, PMF] = floe_dimension(sic);
floe_dim = floe_dim.*1e3;


%- Form drag
tau_f =  0.5*rho_w(:).*(z./floe_dim(:)).*Uice.^2;
%- Friction velocity from form drag
us_f = sqrt(tau_f./rho_w);

% REPLACE NANs in US_IW with 0
us_iw(isnan(us_iw)) = 0;

% Sutherland and Melville (2014) microbreaking coeficient for epsilon
[MB_coef] = Mb_coef_cal(Uwind,sic); %

%- FINAL SHEAR-BASED CALCULATIONS ---%
us_tot = sqrt(Ai.*(us_f.^2 + us_iw.^2) + ((MB_coef'.^(1/3)).*us_aw).^2.*(1-Ai));


%- make a logical index for Obukhov convection criteria
convect = mld(:)./z./Lo;
use = convect<1;
eps = us_tot.^3/.41/z;

%- Calculate Epsilon from heat flux and shear.
%- Reference:
eps_Jb = abs(.58*.84*Jb.*use);
eps_sh = .84*1.76*eps;

%- calculate k from shear (visc. (m2/s), eps (m2/s3)=> m/s
%- multiply by 86400 to get m/d.
k = .419*Sc^(-n).*(vis.*(nansum([eps_Jb,eps_sh],2))).^.25*86400;
keff = (1-Ai).*k;
end



