function [eps_ratio] = Mb_coef_cal(Uwind,f_v)
% function [eps_ratio,wave_age] = Mb_coef_cal(Uwind,f)
% Inputs : Uwind (m/s) , 0<f<1.
for i = 1:length(f_v)
    f = f_v(i);
    x1 = 1E3.*Uwind^2 .* (162.*f.^(-0.49))./9.81; %Smith & Thomson (2016)
    x(i) = x1;
    %%-- Formulation based on Shore protection manual 4th ed by Dept.
    %%-- of the Army, Waterways Experiment Station, Corps of Engineers
    %%--, Coastal Engineering Research Center, Published 1984
    C_D=0.001*(1.1+0.035*Uwind);
    u_f=Uwind.*sqrt(C_D);
    g = 9.8;
    %- fetch limit check
    if 1/(((x1*1E-3).*g/Uwind^2).^(-0.33)*(g./Uwind)) < .398*10^2*u_f/g
        Ts = 1./(((x1*1E-3).*g/Uwind^2).^(-0.33)*(g./Uwind));
    else
        Ts = 2.398*10^2*u_f/g;
    end
    cp = (g/pi/2).*Ts;
    wave_age(i) = cp./Uwind;
    eps_ratio(i) = 21.*(wave_age(i)).^3.5;
end
end