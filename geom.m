function P = geom(P)
%[TR, A, S, cbar, ybar, xbarLE, swp_25, swp_50]
P.e = (1 + P.delta)^(-1); % Span efficiency factor.
P.Cavg = 0.5*(P.cr + P.ct); % Average Chord root
P.TR = P.ct/P.cr; % Taper ratio
P.A = P.b/P.Cavg; % Aspect ratio
P.S = (P.b/2)*P.cr*(1 + P.TR); % Theoretical Wing area
P.cbar = (2/3)*P.cr*((1+P.TR+(P.TR)^2)/(1+P.TR)); % Mean Aerodyanmic chord
P.ybar = (P.b/6)*((1+(2*P.TR))/(1+P.TR)); % Position of Mean Aerodyanmic chord
P.xbarLE = P.ybar * tand(P.swp_LE); % distance between normal wing to the LE of swept wing
P.swp_25 = atand(tand(P.swp_LE) - ((4*(0.25-0)*(1-P.TR))/(P.A*(1+P.TR)))); % Sweep at 25% 
P.swp_50 = atand(tand(P.swp_LE) - ((4*(0.50-0)*(1-P.TR))/(P.A*(1+P.TR)))); % Sweep at 50%
%P.K = P.ao/(2*pi);
%P.a_w = (2*pi*P.A)/(2 + sqrt((((P.A)^2*(P.beta)^2)/(P.K)^2)*(1 + ((tan(P.swp_25/57.3))^2/(P.beta)^2))) + 4); % Lift curve slope

end
