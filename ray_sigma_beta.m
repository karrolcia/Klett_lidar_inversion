% get_ray_sigma_beta

%     wavelength in nm
%     beta and sigma in cm^2/(sr)
%     using data of d.r. bates
%     find the rayleigh extinction and backscattering
%     cross sections (valid from about 300 to 1000 nm)
%     beta and sigma in cm^3/(sr)/km/molecule

%% Atmospheric density calculation
%  av_num = 6.0222e+23;  % molecules/mole
%  rg = 8.314;  % J/mole*K universal gass cnst
%  rho_atm = (P(1:nt,1:nh) .* 1.0e2) ./ (rg .* T(1:nt,1:nh)) .* av_num; % pressure from mb to Pa
%  rho_atm = rho_atm .* 1.0e-6  ; %in 1/cm3

%% 01.05.2014 - this works, magnitued of values is OK!
%%
wavelen = 355 ; % in nm

work = wavelen*1.0e-3 ;
  
 if wavelen < 550
        x = 0.389 * work + 0.09426/work - 0.3228 ;
  else  x = 0.04 ;
 end 
 
sigma = (4.02e-28/ (work ^(4 + x))); % rayleigh extinction coefficient
%backscattering crossection in cm^2sr^-1
% this value is correct according to the D.R. Bates paper


fk = 1.225e-13 * wavelen^4 - ...
     3.911e-10* wavelen^3 + ...
     4.6100e-7 * wavelen^2 - ...
     2.410e-4 * wavelen + 1.095;
 % Effective King Correction factor for Air - calculated correctly

e = (fk - 1.0) * 9.0 / 2.0;

beta_o = 3.0 / (4.0 * pi) * (180.0 + 28.0 * e) / (360.0 + 82 * e) * sigma;

alpha_ray = sigma * 1.0e5; % in cm^3/(sr)/km/molecule
beta_ray = beta_o * 1.0e5 ;% n cm^3/(sr)/km/molecule




