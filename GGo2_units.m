function [o2s]=GGo2_units(T, S,units)

% GGO2 This function calculates the oxygen saturation concentration in umol/kg using the 
%  regression of Garcia and Gordon, 1992,L&O,1307 with Benson and Krause
%  data

% taken from R Hamme's gassat07.m function and adapted to allow calculation
% of o2 sat concentration in umol/kg OR ml/L
% 
% units can equal ml or umol, for ml/L or umol/kg, respectively
% T is in units of C
TS = log((298.15 - T) ./ (273.15 + T));

if ~(strcmp(units,'umol') || strcmp(units, 'ml')),
    error('units term must be either umol or ml')
end

if strcmp(units, 'umol'),
    A0 = 5.80871;
    A1 = 3.20291;
    A2 = 4.17887;
    A3 = 5.10006;
    A4 = -.0986643;
    A5 = 3.80369;
    B0 = -.00701577;
    B1 = -.00770028;
    B2 = -.0113864;
    B3 = -.00951519;
    C0 = -2.75915E-07;
elseif strcmp(units, 'ml'),
    A0 = 2.00907;	
    A1 = 3.22014;
    A2 = 4.05010;
    A3 = 4.94457;
    A4 = -0.256847;	
    A5 = 3.88767;	
    B0 = -0.00624523;	
    B1 = -0.00737614;	
    B2 = -0.010341;
    B3 = -0.00817083;	
    C0 = -0.000000488682;
end



C = A0 + A1 .* TS + A2 .* TS .^ 2 + A3 .* TS .^ 3 + A4 .* TS .^ 4 + ...
    A5 .* TS .^ 5 + S .* (B0 + B1 .* TS + B2 .* TS .^ 2 + B3 .* TS .^ 3) + ...
    C0 .* S .^ 2; 
o2s = exp(C);
end