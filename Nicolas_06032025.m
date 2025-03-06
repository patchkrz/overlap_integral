%% -----
% 
% This code computes field overlap for optical E-field and the E-field
% modulating it. 
%
% * The inputs of this code are optical E-field, modulating E-field,
% electric potential and effective refractive index which is computed at
% the points user defines.
% 
%
% --- Data Export Order ---
% This is important to be able to keep track of data. The following order
% is preferred on COMSOL export interface:
%
% ewfd.Ex [V/m] -> Optical Field x-component
% ewfd.Ez [V/m] -> Optical Field z-component
% es.Ex [V/m]   -> Static Field x-component
% V             -> Electric Potential
% ewfd.neff     -> Effective mode index
%
%
% Results and comparisions:
% https://docs.google.com/spreadsheets/d/11hWJHqptUex4R8OjKD7x7_AgD8aYe10y9XUeD_9BMqI/edit?usp=sharing
%
% -----

data = load('nicolas_06032025_field_amplitudes.txt');

amp_ewfdEx = data(:,23);
amp_ewfdEz = data(:,24);
amp_esEx = data(:,25); 

% NaN values are removed.
validData = ~isnan(amp_ewfdEx) & ~isnan(amp_ewfdEz) & ~isnan(amp_esEx) ;

filtered_ewfdEx = amp_ewfdEx(validData); 
filtered_ewfdEz = amp_ewfdEz(validData);
filtered_esEx = amp_esEx(validData);  

% E = fillmissing(amp_esEx, 'linear');

% meshgrid and multidimensional integration preparetion
[Ex,Ez] = meshgrid(filtered_ewfdEx,filtered_ewfdEz);
ewdfE = sqrt(Ex.^2 + Ez.^2);

n_eff = 1.9286; % effective refractive index (from COMSOL)
v_dc = 1; % control voltage (from COMSOL)
g = 3.2e-6; % gap between signal tracks 
lam = 1055e-9; % operating wavelength
ne = 2.156; % extraordinary refractive index of LN at 1064nm
r33 = 30.8e-12; % [meter/Volt] (Yariv) EO coefficient
L = 4e-3;

nom_overlap = g * trapz(filtered_ewfdEx, trapz(filtered_ewfdEz,(ewdfE.^2).*filtered_esEx,1));

denom_overlap = v_dc * trapz(filtered_ewfdEx,trapz(filtered_ewfdEz,(ewdfE.^2),1));
 
overlap = nom_overlap./denom_overlap;
 
v_pi = (lam*g*n_eff) / (2*ne^4*r33.*overlap);

v_pi_L = v_pi / L;