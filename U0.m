function [Und,etastar,A,B]=U0(Rostar,mustar)
%[Und,etastar,A,B]=U0(Rostar,mustar);
%
% U0:  Given the rossby number (u*/fz) and mustar (u*/fL) calculate the
% nondimensional velocity
%
% Constants Rc=0.2;lamstar=0.028;kappa=0.4
% INPUT Rostar may be a vector
%       mustar may be a vector
% OUTPUT U0 (optional) A B vectors of nondimensional surface velocity
%  actual surface velocity is U0 x u*/etastar
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

% 22 Aug 07 MGM
% Given the rossby number (u*/fz) and mustar (u*/fL) calculate the
% nondimensional velocity
% Constants Rc=0.2;lamstar=0.028;kappa=0.4
% INPUT Rostar may be a vector
%       mustar may be a vector
% OUTPUT U0 (optional) A B vectors of nondimensional surface velocity
%  actual surface velocity is U0 x u*/etastar
lamstar=0.028;
if size(Rostar,2)>size(Rostar,1)
    Rostar=Rostar';
end
if size(mustar,2)>size(mustar,1)
    mustar=mustar';
end
if size(mustar,1)>1 & size(Rostar,1)==1
    Rostar=Rostar*ones(size(mustar));
elseif size(Rostar,1)>1 & size(mustar,1)==1
    mustar=mustar*ones(size(Rostar));
end
Rc=0.2;
kappa=0.4;
delta=(1i/lamstar)^(1/2);
etastar=(1+lamstar*mustar/(kappa*Rc)).^-(1/2);
xisl=-lamstar/kappa;
xi0=(etastar.*Rostar).^-1;
logR=log(Rostar);
logxi=logR+log(etastar*lamstar/kappa);
a=(kappa/lamstar)*(1-etastar)./etastar;
Ue=-1i*delta*exp(delta*xisl);
Und=Ue+(etastar/kappa).*(logxi-(a-delta).*xisl-(delta/2)*a.*xisl.^2);
B=-kappa*imag(Und)./etastar;
A=logR-kappa*real(Und)./etastar;
