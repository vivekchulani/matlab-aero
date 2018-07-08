function Pa = aero(P, Pa, M)

% Calculating Wing Zero Lift angle
[Pa.aol_wo, Pa.aol_m] = Fig_aol_theta(P.swp_25, P.A, P.TR, M, Pa.tovc);
Pa.AOL = Pa.aol + Pa.aol_wo*(P.wo); % Wing zero lift angle
Pa.AOL_m = (Pa.aol + Pa.aol_wo*(P.wo))*Pa.aol_m; % Wing zero lift angle with mach effects included

% Calculating Coeffiecient of Moment of the Aircraft.
[Pa.dCm_theta, Pa.dCm_theta_M] = Fig_dcm_theta(P.swp_25, P.A, P.TR, M);
Pa.Cmac_0 = ((P.A*(cosd(P.swp_25))^2)/(P.A + (2*cosd(P.swp_25))))*Pa.cmac; % Coeffiecient of Moment of the Wing
 
if P.wo == 0
   Pa.CMac = Pa.Cmac_0 * Pa.dCm_theta_M; % Coeffiecient of Moment of the Aircraft without wash out.
else
    Pa.CMac = (Pa.Cmac_0 + (Pa.dCm_theta)*(P.wo))* Pa.dCm_theta_M; % Coeffiecient of Moment of the Aircraft with wash out.
end 

% Calculating aerodynamic centre as a percentage of the chord.
Pa.xac_cr = Fig_xac_cr(P.swp_LE,P.A,P.TR,M);
Pa.xac_cbar = (Pa.xac_cr * P.cr) - P.xbarLE;
Pa.xbar_ac = Pa.xac_cbar/P.cbar;

% Calculating lift curve slope for a finite wing or tail
%Pa.CLa = Pa.ao/(1 + ((Pa.ao*57.3)/(pi*P.e*P.A))); % Lift curve slope Neglecting sweep and mach effects
%Pa.K = (Pa.ao*57.3)/(2*pi);
%Pa.CLa = (2*pi*P.A)/(2 + sqrt((((P.A)^2*(P.beta)^2)/(Pa.K)^2)*(1 + ((tan(P.swp_50/57.3))^2/(P.beta)^2))) + 4); % Lift curve slope
%P.a_t = (2*pi*P.A_t)/(2 + sqrt((((P.A_t)^2*(P.beta)^2)/(P.K)^2)*(1 + ((tan(P.swp_25/57.3))^2/(P.beta)^2))) + 4); % Lift curve slope
%Pa.dwnwashgrad_tail = (2*Pa.a_w)/(pi*Pa.e*Pa.A)
%Pa.CLa = Pa.a_w + Pa.a_t * Pa.nita_t * (S_t/S_w) * (1 - Pa.dwnwash_grad)
Pa.K = (Pa.ao*57.3)/(2*pi);
Pa.CLa = (2*pi*P.A)/(2 + sqrt(((P.A^2*P.beta^2)/Pa.K^2) * (1 + (tand(P.swp_50)^2/P.beta)) + 4)); % Lift curve slope of wing/tail per rad - Subsonic
Pa.CLa = Pa.CLa/57.3;




% Calculating slope of Coeffient of Moment. (Cm_aplha)


end