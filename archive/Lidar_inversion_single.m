% Lidar inversion - module for simulated data
% clc ; clear all ; close all
format long


% load('les_1-merged_with_driz_5m.mat')
% load('true_les_1-merged_with_driz_5m.mat')
load('fire_27_merged_with_mask_v6.mat')
load('fire_27_merged_with_mask_v6_ss.mat')
load('true_fire_27_merged_with_mask_v6.mat')


beta_atten = lid_total_ss ; % Original unit is 1/m/sr, change to 1/km/sr
beta_atten(:,1:2) = 0; % Remove first 3 values - only for the single 
% scattering data - there is an unexplained signal below the cloud
sd_beta_atten = 0.2 .* beta_atten ;

 nh = nz;   
 nt = nx;
 time = linspace(1,nt,nt);
 height = height(1,:) .* 1.0e3;
 height_res = height(2) - height(1);  %Height resolution in km
 

 % Conver Range Corrected Signal to Power
P_sig = zeros(nt,nh);
% P_sig = beta_atten ./ height.^2;

 for it = 1:nt
 P_sig(it,:) = beta_atten(it,:) ./ height.^2;
 end
 sd_P_sig = 0.2 .* P_sig;
 
 para = beta_atten ./ (lin_depol + 1) ;
 perp = beta_atten - para;
 
% Calculate LWP - convert cm to g/m^2
lwp = LWP ./ 1.e-4 ;
lwp = lwp(:,1)';

% Calculate rho & rayleigh scattering
% Constants used in the UV_lidar_inversion
 av_num = 6.0222e+23;  % molecules/mole
 rg = 8.314;  % J/mole*K universal gass cnst
 rho_atm = (P .* 1.0e2) ./ (rg .* T) .* av_num; % pressure from mb to Pa
 rho_atm = rho_atm .* 1.0e-6  ; %in 1/cm3
% Calling ray_sigma_beta script, which provides Rayleigh scattering parameters
 ray_sigma_beta
% Calculating molecular extinction and backscatter fields
 alpha_m = rho_atm .* alpha_ray .*1.0e-3 ;  % in 1/km chnaged to 1/m
 beta_m = rho_atm .* beta_ray  .*1.0e-3 ; % 

  
%% Find peak signal value in each time step

% Find the peak of the attenuated backscatter
 [Peak_P, I_max] = max(P_sig, [], 2) ;
 [Peak_perp, perp_max] = max(perp, [], 2) ;
 
 for icb = 1:nt
%   Peak_P(icb) = max(P_sig(icb,:)) ;
%   I_max(icb)  = find(P_sig(icb,:) >= Peak_P(icb), 1, 'last' );
 
%   Peak_perp(icb) = max(perp(icb,:)) ;    
%   perp_max(icb) = find(perp(icb,:) >= perp(icb), 1, 'last' );
  id_cb_lidar(icb) = find (perp(icb,1:perp_max(icb)) < (Peak_perp(icb)/10), 1,'last');   
  id_cb_radar(icb) = find((Ze(icb,:)) >= 1.e-3, 1, 'first');
  id_ct_radar(icb) = find((Ze(icb,:)) >= 1.e-3, 1, 'last');
 end
 
%% Calculate normalization value for each time step
% Define up and down limit for normalization
down = 6;
up = 11;

lower_bound = I_max+down;
upper_bound = I_max+up;



%% Aerosol Extinction-to-Backscatter ratio used in the inversion process [sr]
% and uncertainity in S used to assess the associated bias in the results 
% extinction-to-backscatter ratio for clouds should be 20
S = zeros(nt,nh);
% S(:,lower_bound:end) = 20.0 ;
% S(:,1:lower_bound-1) = 50.0 ;
S(:,:) = 20.0;
% S = 20.0 ;     
dS = 20.0 ; 
 
for it = 1:nt
      
 %% Normalization height
 % 1. Find the peak value of the signal in each time step

 height_norm_1(it) = height(I_max(it)+down) ;
 height_norm_2(it) = height(I_max(it)+up) ; 
 tau_mol = 0;
 int_S_beta = 0 ; 

%% Calculating Pp, Pp_u, Pp_l and associated standard deviations in a loop
 for iz = 1:nh
 
     
tau_mol = tau_mol + alpha_m(it,iz) * height_res;
int_S_beta = int_S_beta + (beta_m(it,iz) * S(it,iz)) * height_res ;

Pp(iz) = S(it,iz) * P_sig(it,iz) * exp(-2.0*int_S_beta) * exp(2*tau_mol);
   
sd_Pp(iz) = S(it,iz) * sd_P_sig(it,iz) * exp(-2.0 *int_S_beta) * exp(2*tau_mol);

 end
 tau_mmol{it} = tau_mol;
 int_beta{it} = int_S_beta; 
 Pp_int{it} = Pp ;

  
% Calculate the index for io  Find the index corresponding to zo
height_o(it)= (height_norm_1(it) + height_norm_2(it))/2;
io(it)  = find(height <= height_o(it), 1, 'last' );

% Calculating slope of the signal
% Calculating slope of the signal  
x = log(beta_atten(it,I_max(it)+down:I_max(it)+up));
y = height(I_max(it)+down:I_max(it)+up);

beta_slope{it} = x;
height_slope{it} = y;
  
 [poly_coeff , struct_S, mu] = polyfit(beta_slope{it},height_slope{it},1);
 
  slope = poly_coeff(1);

%   [yfit, d_alpha_o] = polyval(poly_coeff,x,struct_S);
%   yresid = y - yfit;
%   SSresid = sum(yresid.^2);
%   SStotal = (length(y)-1) * var(y);
%   rsq = 1 - SSresid./SStotal;

  
  lm = fitlm(x,y);
alpha_o_p(it,:) = -0.5 * (diff(beta_slope{it}) ./ height_res) +  beta_m(it,io(it)) * S(it,io(it));
alpha_slope(it) = alpha_o_p(it,length(beta_slope{it})/2); 
% alpha_slope(it) = alpha_o_p(it,1); 
% The io is in the middle of the normalisation interval

  
% Estimating value of alpha_o 

alpha_o(it) = Ext(it,io(it)) + S(it, io(it)) .* beta_m(it, io(it));
% alpha_o_p(it) = 0.2;
%  alpha_o_p(it) = -0.5 * (slope(it));% + (1/(Pp_int{it}(io(it)))+ beta_m(it,io(it)) * S(it,iz)));

[alphap{it}, sd_alphap{it}] = Klett_inversion(height, Pp_int{it}, sd_Pp, ...
     io(it), alpha_slope(it)) ;

alpha(it,:) = alphap{it} ;%- S(it,:) .* beta_m(it,:) ; % becuse there is no aerosols below the cloud in ECSIM simulations
sd_alpha(it,:) = sd_alphap{it} ;%- S(it,:) .* beta_m(it,:) ;

   
 % Calculating the corresponding beta (backscattering) values
beta(it,:) = alpha(it,:) ./ S(it,iz) ;
sd_beta(it,:) = sd_alpha(it,:) ./ S(it,iz) ;
  
clear Pp sd_Pp x y
end
       

tau = height_res * cumsum(alpha,1);         %Optical depth
sd_tau = height_res * cumsum(sd_alpha,1) ;  %Estimated random error in optical depth
e = std(alpha);

tt=50;

figure('NumberTitle','off', ...
    'Units', 'centimeters','Position',[2 20 15 15])
% errorbar(alpha(100,:),e)
plot(alpha(tt,1:30), height(1:30), '.-r')
% herrorbar(alpha(nt,:),sd_alpha(nt,:))
hold on
plot(Ext(tt,1:30), height(1:30), '.-b')
hold on 
normalisation = patch([0.25 -0.1 -0.1 0.25], ...
    [height(I_max(tt)+down) height(I_max(tt)+down) height(I_max(tt)+up) height(I_max(tt)+up)], ...
    'b');
set(normalisation,'FaceColor','b', ...
    'FaceAlpha', 0.1)
% refline([0 height(I_max(tt))])
refline([0 height(id_cb_radar(tt))])
hold on
% plot(alpha_corr(1,:), height, '*-r')
%refline(0,height(1,ice),'r')
% cloud_base = refline(0,height(1,icb));
% set(cloud_base,'Color','red')
% title('Retrieved Extinction coefficient');
 