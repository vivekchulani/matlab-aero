function xac_cr = Fig_xac_cr(SweepLE,A,Taper,M)
%xac_cr = Fig_xac_cr(Sweep_LE, Aspect Ratio, Taper Ratio, Mach Number)
%Aerodynamic Center with respect to chord root, Datcom Figure 4.1.4.2-26
%Input
%SweepLE = Leading edge sweep
%Aspect Ratio
%TR = Taper ratio
%M = Mach
%Output
%xac_cr = Aerodynamic center location on the chord root in chord roots

%Taper Ratio
X3 = [0., .2, .25, .33, .5, 1.];
%Aspect Ratio
X2 = [1., 2., 3., 4., 5., 6.];
%c/4 Sweep
X1 = [0., .2, .4, .6, .8, 1.];
%Subsonic Range 0 - 1
Y1 = {{[ .250, .245, .240, .235, .230, .225] ...
       [ .335, .335, .335, .335, .335, .335] ...
       [ .420, .430, .435, .445, .450, .455] ...
       [ .500, .515, .530, .540, .550, .560] ...
       [ .580, .600, .630, .645, .660, .680] ...
       [ .680, .695, .720, .740, .760, .780]};
 
      {[ .285, .275, .270, .265, .260, .255] ...
       [ .400, .410, .415, .415, .415, .410] ...
       [ .510, .530, .535, .540, .545, .550] ...
       [ .640, .650, .660, .675, .685, .690] ...
       [ .750, .765, .780, .785, .800, .815] ...
       [ .870, .880, .895, .905, .920, .930]};

      {[ .300, .295, .285, .280, .275, .265] ...
       [ .420, .420, .425, .425, .425, .430] ...
       [ .545, .550, .560, .565, .575, .580] ...
       [ .670, .680, .690, .700, .710, .720] ...
       [ .795, .805, .815, .830, .840, .850] ...
       [ .925, .945, .960, .965, .975, .980]};
 
      {[ .325, .320, .315, .305, .300, .290] ...
       [ .460, .460, .460, .460, .455, .455] ...
       [ .595, .600, .600, .600, .610, .620] ...
       [ .735, .740, .750, .760, .765, .775] ...
       [ .885, .890, .895, .900, .910, .925] ...
       [1.045,1.050,1.050,1.060,1.065,1.075]};
 
      {[ .355, .350, .345, .340, .330, .320] ...
       [ .530, .530, .525, .525, .520, .520] ...
       [ .700, .700, .700, .705, .710, .710] ...
       [ .880, .880, .885, .890, .890, .895] ...
       [1.040,1.045,1.050,1.055,1.060,1.065] ...
       [1.200,1.205,1.210,1.215,1.225,1.230]};
       
      {[ .510, .490, .480, .470, .460, .450] ...
       [ .750, .750, .750, .750, .745, .740] ...
       [1.000,1.000,1.000,1.000,1.000,1.000] ...
       [1.250,1.250,1.250,1.250,1.250,1.250] ...
       [1.500,1.500,1.500,1.490,1.490,1.490] ...
       [1.740,1.740,1.740,1.730,1.730,1.730]}};

% Subsonic Range 1 - 0 
Y2 = {{[ .165, .180, .200, .210, .220, .225] ...
       [ .335, .335, .335, .335, .335, .335] ...
       [ .500, .480, .465, .460, .460, .455] ...
       [ .670, .625, .595, .580, .575, .560] ...
       [ .830, .750, .730, .705, .695, .680] ...
       [ .990, .860, .835, .810, .795, .780]};

      {[ .200, .215, .230, .240, .250, .255] ...
       [ .400, .400, .400, .405, .410, .410] ...
       [ .600, .580, .565, .560, .555, .550] ...
       [ .795, .760, .735, .715, .700, .690] ...
       [ .970, .910, .870, .840, .825, .815] ...
       [1.150,1.050,1.000, .965, .940, .930]};
 
      {[ .230, .240, .245, .250, .260, .265] ...
       [ .415, .420, .425, .425, .430, .430] ...
       [ .630, .615, .600, .590, .585, .580] ...
       [ .830, .785, .760, .740, .730, .720] ...
       [1.030, .950, .905, .880, .865, .850] ...
       [1.250,1.090,1.050,1.015, .990, .980]};
 
      {[ .220, .240, .250, .265, .280, .290] ...
       [ .440, .445, .450, .450, .455, .455] ...
       [ .670, .655, .640, .630, .625, .620] ...
       [ .880, .830, .805, .790, .780, .775] ...
       [1.070,1.000, .960, .940, .935, .925] ...
       [1.270,1.170,1.120,1.100,1.085,1.075]};
 
      {[ .250, .270, .295, .310, .315, .320] ...
       [ .500, .505, .510, .515, .520, .520] ...
       [ .750, .740, .730, .720, .715, .710] ...
       [ .980, .940, .915, .900, .900, .895] ...
       [1.190,1.120,1.090,1.080,1.070,1.065] ...
       [1.380,1.300,1.270,1.250,1.240,1.230]};
 
      {[ .340, .380, .410, .430, .440, .450] ...
       [ .680, .700, .720, .730, .730, .740] ...
       [ .950, .980,1.000,1.000,1.000,1.000] ...
       [1.200,1.230,1.250,1.250,1.250,1.250] ...
       [1.440,1.470,1.480,1.480,1.490,1.490] ...
       [1.680,1.710,1.710,1.720,1.720,1.730]}};

%Supersonic Range 0 - 1
Y3 = {{[ .165, .210, .250, .290, .310, .345] ...
       [ .335, .365, .390, .415, .445, .470] ...
       [ .500, .540, .560, .560, .560, .560] ...
       [ .670, .670, .670, .670, .670, .670] ...
       [ .830, .775, .775, .775, .775, .775] ...
       [ .990, .930, .895, .895, .895, .895]};
 
      {[ .200, .230, .280, .305, .335, .360] ...
       [ .400, .445, .485, .500, .520, .530] ...
       [ .600, .630, .650, .660, .665, .665] ...
       [ .795, .800, .800, .805, .810, .815] ...
       [ .970, .965, .955, .955, .955, .955] ...
       [1.150,1.135,1.120,1.100,1.100,1.105]};
 
      {[ .230, .275, .300, .330, .350, .370] ...
       [ .415, .470, .500, .530, .545, .550] ...
       [ .630, .670, .680, .685, .690, .690] ...
       [ .830, .835, .835, .840, .845, .850] ...
       [1.030,1.015,1.005,1.000,1.005,1.010] ...
       [1.250,1.225,1.200,1.170,1.165,1.160]};

      {[ .220, .280, .315, .345, .375, .390] ...
       [ .440, .500, .535, .560, .570, .580] ...
       [ .670, .700, .720, .725, .740, .740] ...
       [ .880, .885, .895, .900, .900, .900] ...
       [1.070,1.070,1.075,1.075,1.080,1.080] ...
       [1.270,1.260,1.260,1.255,1.255,1.255]};

      {[ .250, .300, .330, .380, .415, .445] ...
       [ .500, .560, .600, .620, .635, .640] ...
       [ .750, .780, .800, .820, .820, .825] ...
       [ .980, .990,1.000,1.020,1.020,1.020] ...
       [1.190,1.200,1.200,1.210,1.220,1.225] ...
       [1.380,1.390,1.400,1.410,1.420,1.420]};
 
      {[ .340, .380, .410, .460, .500, .540] ...
       [ .680, .700, .730, .770, .790, .840] ...
       [ .950, .990,1.010,1.050,1.080,1.120] ...
       [1.200,1.240,1.290,1.330,1.370,1.420] ...
       [1.440,1.500,1.550,1.610,1.670,1.720] ...
       [1.680,1.760,1.820,1.890,1.950,2.020]}};

%Supersonic Range 1 - 0
Y4 = {{[ .415, .410, .400, .385, .370, .345] ...
       [ .500, .500, .495, .485, .480, .470] ...
       [ .585, .580, .580, .575, .570, .560] ...
       [ .670, .670, .670, .670, .670, .670] ...
       [ .750, .750, .755, .760, .765, .775] ...
       [ .830, .840, .845, .855, .870, .895]};

      {[ .460, .455, .445, .420, .390, .360] ...
       [ .575, .575, .570, .560, .545, .530] ...
       [ .695, .695, .690, .685, .680, .665] ...
       [ .800, .805, .805, .810, .815, .815] ...
       [ .920, .930, .935, .945, .970, .955] ...
       [1.040,1.045,1.050,1.075,1.110,1.105]};

      {[ .475, .465, .450, .430, .400, .370] ...
       [ .600, .600, .595, .585, .575, .550] ...
       [ .725, .730, .730, .725, .715, .690] ...
       [ .850, .850, .855, .865, .870, .850] ...
       [ .970, .975, .980,1.000,1.020,1.010] ...
       [1.110,1.110,1.110,1.130,1.180,1.160]}; 

      {[ .500, .490, .470, .450, .425, .390] ...
       [ .640, .635, .630, .620, .600, .580] ...
       [ .770, .775, .780, .775, .765, .740] ...
       [ .920, .915, .920, .930, .935, .900] ...
       [1.050,1.055,1.060,1.080,1.105,1.080] ...
       [1.195,1.200,1.205,1.225,1.265,1.255]};

      {[ .550, .535, .525, .500, .475, .445] ...
       [ .720, .715, .710, .690, .670, .640] ...
       [ .890, .890, .890, .885, .870, .825] ...
       [1.060,1.050,1.050,1.060,1.065,1.020] ...
       [1.215,1.215,1.220,1.245,1.270,1.225] ...
       [1.380,1.380,1.395,1.420,1.470,1.420]};

      {[ .760, .730, .700, .650, .600, .540] ...
       [1.000,1.000, .970, .930, .890, .840] ...
       [1.240,1.230,1.230,1.220,1.190,1.120] ...
       [1.500,1.480,1.480,1.490,1.470,1.420] ...
       [1.750,1.720,1.730,1.760,1.780,1.720] ...
       [2.000,1.970,1.980,2.020,2.070,2.020]}};


nmY = length(X3);
nmX2 = length(X2);
Z1 = cell(size(nmY));
Z2 = zeros(size(nmX2));
AT = A*tand(SweepLE);
if M <= 1;
   B = sqrt(1 - M^2);
   if B == 0;
      X0 = 1; Y0 = Y2;
   else
      X0 = tand(SweepLE)/B;
      if X0 > 0 || X0 <= 1; Y0 = Y1; end
      if X0 > 1; X0 = 1/X0; Y0 = Y2; end
   end
end
if M >= 1;
   B = sqrt(M^2 - 1);
   if B == 0;
      X0 = 1; Y0 = Y3;
   else
      X0 = tand(SweepLE)/B;
      if X0 > 0 || X0 <= 1; Y0 = Y4; end
      if X0 > 1; X0 = 1/X0; Y0 = Y3; end
   end
end
for i=1:nmY
    for j=1:nmX2
        %interpolate by sweep for each graph under each aspect ratio
        %and taper ratio
        Z1{i}(j) = interp1(X1,Y0{i}{j},X0,'spline','extrap');
    end
    %interpolate by aspect ratio for each graph under taper ratio
    Z2(i) = interp1(X2,Z1{i},AT,'linear','extrap');
end
%interpolate by taper ratio
xac_cr = interp1(X3,Z2,Taper,'spline');

end

