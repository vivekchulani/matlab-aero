clear;
clc;

FC.vel = 175; %input('\n At what Velocity is the Aircraft flying?: ');
FC.alt = 8000; %input('\n At what Altitude is the Aircraft flying?: ');
FC.wt = 1700; %input('\n What is the Weight of the Aircraft?: ');

% Calculating Atmospheric Properties
if (FC.alt <= 278386) && (FC.alt >= 0)
    
    [FC.T, FC.P, FC.De, FC.Sos, FC.Dv, FC.Kv] = atmos(FC.alt);
    %fprintf('\n The Standard Atmospheric conditions at %.2f ft are as below: \n Temperature: %.4f °R \n Pressure: %.4f lbf/ft2 \n Density: %.4f slug/ft3\n Speed of Sound: %.4f ft/s \n Dynamic Viscosity: %.4f slug/(s·ft) \n', FC.alt, FC.T, FC.P, FC.De, FC.Sos, FC.Dv);
    
else
    
    msgbox('This Program cannot determine standard atmospeheric conditions this high.')
    
end

FC.M = FC.vel/FC.Sos; % Calculating Mach Number

if FC.M <= 0.8 % Calulating Compressibility effects designated by beta
    WG.beta = sqrt(1 - (FC.M)^2);
else
    WG.beta = sqrt((FC.M)^2 - 1);
end




% Wing Section
WG.delta = 0.07; 
WG.b = 33;         % Wing Span in feet
WG.cr = 5;  % Chord Root in feet
WG.cr_e = 5; % Chord Root exposed in feet
WG.ct = 5;   % Chord Tip in feet
WG.swp_LE = 0;  % Sweep Angle at LE in degrees
WG.wo = 0;  % Washout in degrees
WG.i = 0; % Wing incidence in degrees
WG.XW = 6.25; % Nose to LE of wing

WG  = geom(WG);

% Calculating Aerodynamic Parameters (Airfoil parameters) of the Wing
%WGa.aol_cr = -2; % Zero lift angle of the airfoil at chord root
%WGa.aol_ct = 0; % Zero lift angle of the airfoil at chord tip
%WGa.aol = 0.5*(WGa.aol_cr + WGa.aol_ct); % Average Zero lift angle
WGa.aol = -2.8;
%WGa.ao_root = 0.104; % 2-D lift curve slope for airfoil at the root
%WGa.ao_tip = 0.114; % 2-D lift curve slope for airfoil at the tip
%WGa.ao = 0.5*(WGa.ao_root + WGa.ao_tip); % 2-D lift curve slope for combined airfoils 
WGa.ao = .118;
%WGa.cmac_root = -0.047; % Coefficient of Moment about the mean aerodynamic centre - root
%WGa.cmac_tip = 0; % Coefficient of Moment about the mean aerodynamic centre - tip
%WGa.cmac = 0.5*(WGa.cmac_root + WGa.cmac_tip); % Coefficient of Moment about the mean aerodynamic centre - combined
WGa.cmac = -0.069;
%WGa.ac_root = 0.247; % Aerodynamic centre - root
%WGa.ac_tip = 0.25; % Aerodynamic centre - tip
%WGa.ac =  0.5*(WGa.ac_root + WGa.ac_tip); % Aerodynamic centre - combined
WGa.ac = .262;
%WGa.mcr_root = 0.69; % Critical Mach number - root
%WGa.mcr3D_root = WGa.mcr_root/cosd(WG.swp_LE); % Critical Mach number 3D - root
%WGa.mcr_tip = 0.727; % Critical Mach number - tip
%WGa.mcr3D_tip = WGa.mcr_tip/cosd(WG.swp_LE); % Critical Mach number 3D - tip
%WGa.mcr = 0.5*(WGa.mcr_root + WGa.mcr_tip); % Critical Mach number - combined
%WGa.mcr3D = 0.5*(WGa.mcr3D_root + WGa.mcr3D_tip);
%WGa.mcr3D = 0;
%WGa.tovc_root = 0.12; %*WG.cr; % Thickness chord ratio
%WGa.tovc_tip =  0.12; %*WG.ct; % Thickness chord ratio
%WGa.tovc = 0.5*(WGa.tovc_root + WGa.tovc_tip); % Thickness chord ratio
WGa.tovc = 0.15;

WGa = aero(WG, WGa, FC.M);
%WGa.K = (WGa.ao*57.3)/(2*pi);
%WGa.CLa = (2*pi*WG.A)/(2 + sqrt(((WG.A^2*WG.beta^2)/WGa.K^2) * (1 + (tan(WG.swp_50/57.3)^2/WG.beta)) + 4)); % Lift curve slope of wing per rad
%WGa.CLa = WGa.CLa/57.3;
%WG.dwnwashgrad = (WGa.CLa)/(pi*WG.e*WG.A);


% Horizontal tail section
HT.delta  = 0.01;
HT.b = 8; % Span in feet
HT.cr = 4; % Chord Root in feet
HT.ct = 4; % Chord Tip in feet
%HT.A = 4.85; % Aspect Ratio
%HT.a_t = 0.077; % Lift curve slope
HT.swp_LE = 0; % Sweep at LE
%HT.St = 32.2; % Area
HT.wo = 0;
HT.beta = WG.beta;
HT.eita = 0.9;
HT.i = -1.5; % Horizontal Tail incidence


HT = geom(HT);
%HT.dwnwashgrad = 2*WG.dwnwashgrad;
% Subsonic
HT.l_h = 14.3; %input(); % wing mac c/4 to tail mac c/4 parallel to wing root chord
HT.h_h =  4.2; %input(); % Height of tail mac c/4 above (+) or below (-) the wing mac c/4 perpendicular to wing root chord
HT.K_A = (1/WG.A) - (1/(1 + WG.A^2)); % Correction factor
HT.K_TR = (10 - 3*WG.TR)/7; % Correction factor
HT.K_H = (1 - (HT.h_h/WG.b))/(2*(HT.l_h/WG.b))^(1/3); % Correction factor
HT.dwnwashgrad = 4.44*(HT.K_A * HT.K_TR * HT.K_H*(cosd(WG.swp_25))^0.5)^1.19; % Downwash Gradient at horizontal tail - Subsonic
%HT.dwnwashgrad = 0.1;
% Supersonic
%HT.dwnwashangle = (1.62*WGa.CL)/(pi*WG.A*WG.e);
%HT.dwnwashgrad = (1.62*WGa.CLa)/(pi*WG.A*WG.e);

% Calculating Aerodynamic Parameters (Airfoil parameters) of the Horizontal tail
HTa.aol = -2.8; % Zero lift angle
HTa.ao = 0.118; % 2-D lift curve slope
HTa.cmac = -0.069; % Coefficient of Moment about the mean aerodynamic centre
HTa.ac = 0.262; % Aerodynamic centre
HTa.mcr = 0.727; % Critical Mach number
HTa.mcr3D_root = HTa.mcr/cosd(HT.swp_LE);
HTa.tovc = 0.15;

HTa = aero(HT, HTa, FC.M);



% Aircraft
AC.XCG = 7.6; % Distance from nose to CG location
AC.l_t = 14.0; % distance from cg to a.c. of HT
AC.l_f = WG.XW + 0.25*WG.cr; % wing's 1/4 root chord position on fuselage
AC.L_f = 11.5; % Length of fuselage
AC.xbar_CG = 0.27; % Non dimensionalised with WG.cr
AC.W_fmax = 4; % Maximum width of fuselage
%AC.W_nmax = 2.8; % Mean width at each segment
%AC.xbar_n = 3.6; % Length of segment forward of chord root exposed.
%AC.dEu_da = -1 + 1.69*(AC.xbar_n/WG.cr_e)^(-0.4543); % Variation of upwash with distance from wing leading edge
%AC.dEu_da_c = AC.dEu_da*(WGa.CLa/0.0785); % Corrected 
%AC.db_da = ( 1 + AC.dEu_da_c); % Local upwash the fuselage experiences as a function of alpha

% Lift Curve slope of the AC
AC.CLapha = WGa.CLa + HTa.CLa * HT.eita * (HT.S/HT.S) * (1 - HT.dwnwashgrad);

% Fuselage
AC.K_f = 0.002322*exp(5.037*(AC.l_f/AC.L_f)); % Emperical factor in per degree
AC.CM_alpahaf = (AC.K_f * (AC.W_fmax)^2 * AC.L_f)/(WG.S * WG.cbar); % Gilruth & White method
AC.CM_CLf = AC.CM_alpahaf/WGa.CLa; % Stability of fuselage

% Wing
AC.CM_alpahaw = WGa.CLa*(AC.xbar_CG - WGa.xbar_ac);

% Tail
AC.Vbar_t = (HT.S * AC.l_t)/(WG.S * WG.cbar);
AC.CM_alpahat = -(HTa.CLa * HT.eita * AC.Vbar_t * (1 - HT.dwnwashgrad));

% Nacelles
AC.CM_alpahan = 0;
%AC.CM_alpahan = ((pi/2)*((AC.W_nmax^2*AC.db_da*AC.xbar_n)/(WG.S*WG.cbar)))*2*(1/57.3); % Per degree

% CM_alpha AC
AC.CM_alpahaAC = AC.CM_alpahaf + AC.CM_alpahaw + AC.CM_alpahat + AC.CM_alpahan;


