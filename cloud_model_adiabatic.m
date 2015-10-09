% Adiabatic formulas for cloud model
% Based on Donovan et al., 2014

% Parameters that need to be set include: 
% - cloud base height
% - Effective radius at 100 m in the cloud [micrometers]
% - LWC' [gm^-3 * km-1] - vertical derivative of the LWC

format long

%% Load height and density profile and create time series. 

load('rho_cld.mat') % in cm^-3
rho_atm = rho_cld*1.0e6; 
% in m^-3 this should be value in the range 1.0e25
rho = 1000; % density of water in kg/m^3
load('height.mat') % in meters 
height = height';
% Transpose matrix to have same dimensions for all the variables
% define time & height

time = linspace(0,2,120);

nt = length(time);
nz = length(height);

%% Definie necessary parameters
%r_eff_100 = 2.0 ; % Effective radius at 100 m in the cloud [in micrometers]
% A time series of effective radius (values from 2 to 5 micrometers)
r_eff_100 = linspace(4.0e-6,2.0e-6,nt) ; % in meters
% Cloud height
icb = 70; % index where the cloud begins
ice = 100; % index where the cloud finishes
cb = height(icb) ; % cloud base height
LWC_der = 0.2 * 1.0e-6 ; %  LWC' [kg / m^3 * m] - vertical derivative of the LWC


%% 1. Calculate adiabatic data profiles 
%Cloud droplet effective radius; (icb+7) - effective radius 100m from the 
%cloud base.

% 1. Cloud droplet effective radius;
r_eff = zeros (nt,nz);
% 2. Extinction coefficient profile
Alpha = zeros(nt,nz) ;

for it = 1:nt
    % 1. Cloud droplet effective radius;
    r_eff(it,icb:ice) = (r_eff_100(it) / 100 ^(1/3)) .* (height(icb:ice) - height(icb)).^(1/3) ;
    % 2. Extinction coefficient profile
    Alpha(it,icb:ice) = 3/2 * (100^(1/3) / rho) * (LWC_der/ r_eff_100(it)) * (height(icb:ice) - height(icb)).^(2/3);
end

% for it = 1:nt
%     for j = icb:ice
%     % 1. Cloud droplet effective radius;
%     r_eff(it,j) = (r_eff_100(it) / 100 ^(1/3)) * (height(j) - height(icb))^(1/3) ;
%     % 2. Extinction coefficient profile
%     Alpha(it,j) = 3/2 * (100^(1/3) / rho) * (LWC_der/ r_eff_100(it)) * (height(j) - height(icb))^(2/3);
% 
%     end
% end

%% Calculate the molecular extinction and backscatter coefficients
% Calling ray_sigma_beta script, which provides Rayleigh scattering parameters
ray_sigma_beta
% Calculating molecular extinction and backscatter fields
alpha_m = rho_atm * alpha_ray ;  % now in 1/m - those values are correct (should be 1.0e-10)
beta_m = rho_atm * beta_ray   ; % now in 1/m
z_res = height(2) - height(1);

tau_mol = zeros(nt,nz);
tau_cld = zeros(nt,nz);
for it = 1:nt
    % Calculate molecular optical thickness
    tau_mol(it,:)= z_res * cumsum(alpha_m(it,:));
    % Calculate cloud optical thickness
    tau_cld(it,:) = z_res * cumsum(Alpha(it,:));
end

alpha_total = alpha_m ;
alpha_total(:,icb:ice) = Alpha(:,icb:ice) ; 


%% Calculate attenuated backscatter
% Aerosol Extinction-to-Backscatter ratio used in the inversion process [sr]
% and uncertainity in S used to assess the associated bias in the results 
% extinction-to-backscatter ratio for clouds should be 20
S = 20.0 ; 
dS = 0.1 * S ; 

beta_atten = (alpha_total./S + beta_m) .* exp(-2.0 .* (tau_cld + tau_mol)) ;

sd_beta_atten = 0.1 * beta_atten ;

%beta_noise = addAWGN(beta_atten, noise);

beta_noise = awgn(beta_atten,30,'measured');

% for it = 1:nt,
% 
%   % Change from attenuated backscatter to power  
%   P(it,:) = beta_atten(it,:) .* height.^2 ;
%   sd_P(it,:) = sd_beta_atten(it,:) .* height.^2 ;
%   
% end

%% Plot results - Attenuaed backscatter time series depended on the height
figure
plot (beta_atten, height)
title('Attenuated backscatter');
figure
plot (alpha_total, height)
title('Modelled Extinction coefficient');
plot (beta_noise, height)
title('Attenuated backscatter with noise');
% figure
% plot (P, height)
% title('Range Corrected Signal');
