function [T, P, De, Sos, Dv, Kv] = atmos(alt)
%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%                                                                        %
% Function calculates the standard atmospheric properties                %
% In the first 2 layers of the atmoshere.                                %
% Latim Jonathan ERAU                                                    %
%                                                                        %
% Symbols: T: Temperature (°R)                                           %
%          P: Pressure (lbf/ft^2)                                        %
%          De: Density (slug/ft^3)                                       %
%          Sos: Speed of Sound (ft/s)                                    %
%          Dv: Dynamic Viscosity (slug/(s·ft))                           %
%          Kv: Kinematic Viscosity (ft^2/s)                              %
%                                                                        %
% Input: Altitude at which one desires to calculate standard             %
%        Atmospheric properties in the first 2 layers of the atmoshere.  %
% Output: Temperature, Pressure, Density, Speed of Sound, Dynamic        %
%         Viscosity, Kinematic Viscosity.                                %
%                                                                        %
%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

t =  518.69; % Temperature at sea level
p = 2116 ; % Pressure at sea level
de = 0.002377; % Density at sea level
s = 1116; % Speed of sound at sea level
dv = 3.737*10^-7; % Dynamic Viscosity at sea level
kv = 1.46072 * 10^(-5); % Kinematic Viscosity at sea level

if alt <= 36089
    T = t*(1-(alt/145442));
    P = p*((1-(alt/145442)))^5.255876;
    De = de*((1-(alt/145442)))^4.255876;
    
    
elseif (alt <= 65617) && (alt > 36089) % Isothermal
    T = t*0.751865;
    P = p*0.223361*exp(-(alt-36089)/20806);
    De = de*0.297076*exp(-(alt-36089)/20806);
    
    
elseif (alt <= 104987) && (alt > 65617) % Inversion
    T = t*(0.682457 + (alt/945374));
    P = p*(0.988626 + (alt/652600))^(-34.16320);
    De = de*(0.978261 + (alt/659515))^(-35.16320);
    
    
elseif (alt <= 154199) && (alt > 104987) % Inversion
    T = t*(0.482561 + (alt/337634));
    P = p*(0.898309 + (alt/181373))^(-12.20114);
    De = de*(0.857003 + (alt/190115))^(-13.20114);
    
    
elseif (alt <= 167323) && (alt > 154199) %Isothermal
    T = t*(0.939268);
    P = p*(0.00109456*exp(-(alt-154199)/25992));
    De = de*(0.00116533*exp(-(alt-154199)/25992));
    
    
elseif (alt <= 232940) && (alt > 167323) % see ERRATA
    T = t*(1.434843 - (alt/337634));
    P = p*(0.838263 - (alt/577922))^(12.20114);
    De = de*(0.798990 - (alt/606330))^(11.20114);
    
    
elseif (alt <= 278386) && (alt > 232940) % see ERRATA
    T = t*(1.237723 - (alt/472687));
    P = p*(0.917131 -(alt/637919))^(17.08160);
    De = de*(0.900194 -(alt/649922))^(16.08160);
    
    
end
Sos = s*(T/t)^0.5;
Dv = dv*(T/t)^0.7;
Kv = kv*(T/t)^(-3.25892);

end
