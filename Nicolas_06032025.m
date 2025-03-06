data = load('nicolas_06032025_field_amplitudes.txt');

amp_ewfdEx = data(:,23);
amp_ewfdEz = data(:,24);
amp_esEx = data(:,25); 

% NaN values are removed.
validData = ~isnan(amp_ewfdEx) & ~isnan(amp_ewfdEz) & ~isnan(amp_esEx) ;

filtered_ewfdEx = amp_ewfdEx(validData); 
filtered_ewfdEz = amp_ewfdEz(validData);
filtered_esEx = amp_esEx(validData);  

% --- Coordinate values ---
% It is the crossection area of the system. Simulator calls dimensions as x
% and y-coordinates.    
% x = data(:,1);
% y = data(:,2);
% 
% x = x(validData); 
% y = y(validData); 

% meshgrid and multidimensional integration preparetion
[Ex,Ez] = meshgrid(filtered_ewfdEx,filtered_ewfdEz);
ewdfE = sqrt(Ex.^2 + Ez.^2);

n_eff = 1.9286; % effective refractive index
v_dc = 1; % control voltage
g = 3.2e-6; % gap between signal tracks 
lam = 1055e-9; % operating wavelength
ne = 2.156; % extraordinary refractive index of LN at 1064nm
r33 = 30.8e-12; % [meter/Volt] (Yariv) EO coefficient
L = 1e-2;

nom_overlap = g * trapz(filtered_ewfdEx, trapz(filtered_ewfdEz,(ewdfE.^2).*filtered_esEx,1));
% deneme = trapz(filtered_ewfdEx,(ewdfE.^2).*filtered_esEx,1);
% deneme_devam = g * trapz(filtered_ewfdEz, deneme);

denom_overlap = v_dc * trapz(filtered_ewfdEx,trapz(filtered_ewfdEz,(ewdfE.^2),1));
 
overlap = nom_overlap./denom_overlap;
 
v_pi = (lam*g*n_eff) / (2*ne^4*r33.*overlap);

v_pi_L = v_pi / L;