% Lidar inversion - module for simulated data
 %% Read in the adiabatic data
clc ; clear all ; close all

% load('fire_27_merged_with_mask_v6.mat')
% load('true_fire_27_merged_with_mask_v6.mat')
% load('fire_27_merged_with_mask_v6_ss.mat')

 load ('les_1-merged_with_driz.mat' )
 load('true_les_1-merged_with_driz.mat')
 
%  cloud_model_adiabatic
 beta_atten = lid_total ;% Use he signal with noise
 
 % For the "fire" simulation. In the "les" simulations those values should
% be reverted
 nh = nz;
 nt = nx;
 time = linspace(1,nt,nt);
 height = height(1,:) .* 1.e3;
 height_res = height(2) - height(1);  %Height resolution in km


 %perp = (beta_atten .* lin_depol) ./ 1 + lin_depol ;
 para = beta_atten ./ (lin_depol + 1) ;
 perp = beta_atten - para;
%   beta_atten2 = para + perp;

 sd_beta_atten = 0.2 .* beta_atten ;

% Choose vartical extent of interest. Define max altitude in meters 

alt_min = 0.0 ;
alt_max = 1000.0 ;
alt = height(height<alt_max);
% Choose data only in the vertical extent of interest (except for the
% beta_atten - that should be done after the inversion is complete.

% Ze (Ze == 0 & Ze == -0) = NaN;
% Ze = Ze(:, height<alt_max);
% Vd = Vd (:, height<alt_max);
% LWP = LWP (:, height<alt_max);
% No = No (:, height<alt_max);
% Ext = Ext (:, height<alt_max);
% LWC = LWC (:, height<alt_max);

% Calculate LWP - convert cm to g/m^2
lwp = LWP ./ 1.e-4 ;
lwp = lwp(:,1)';

% Calculate rho & rayleigh scattering
% Constants used in the UV_lidar_inversion
 av_num = 6.0222e+23;  % molecules/mole
 rg = 8.314;  % J/mole*K universal gass cnst
 rho_atm = (P .* 1.0e2) ./ (rg .* T) .* av_num; % pressure from mb to Pa
% Calling ray_sigma_beta script, which provides Rayleigh scattering parameters
 ray_sigma_beta
% Calculating molecular extinction and backscatter fields
 alpha_m = rho_atm .* alpha_ray ;  % now in 1/m - those values are correct (should be 1.0e-10)
 beta_m = rho_atm .* beta_ray   ; % now in 1/m
 z_res = height(2) - height(1);
 tau_mol = z_res .* cumsum(alpha_m);

  
%% Find peak signal value in each time step

% Find the peak of the perpendicular signal 
 [Peak_perp, perp_max] = max(perp, [], 2) ;
% Find the peak of the attenuated backscatter
 [Peak_P, I_max] = max(beta_atten, [], 2) ;
 
% Find the height where the signal is 1/10 of the Peak perpendicular value
  
 for icb = 1:nt
  id_cb_lidar(icb) = find (perp(icb,1:perp_max(icb)) < (Peak_perp(icb)/10), 1,'last');
%   id_ct_lidar(icb) = find (perp(icb,perp_max(icb):I_max(icb)) ~= 0, 1,'last');
  id_cb_radar(icb) = find((Ze(icb,:)) >= 1.e-3, 1, 'first');
  id_ct_radar(icb) = find((Ze(icb,:)) >= 1.e-3, 1, 'last');
 end
 
%% Calculate normalization value for each time step
% Calculate slope of the regression line - use the value as alpha_o
% slope = zeros(1,nt);         % Slope of the regression line
% slope_2 = zeros(1,nt);       % Slope of the regression line - 2nd degree polynomial
% Define up and down limit for normalization
down = 5;
up = 8;

%% Corecting signal for the multiple scattering
perp_int = zeros(nt,nh);
para_int = zeros(nt,nh);

lower_bound = id_cb_lidar;
upper_bound = I_max+up+6;
perp_sum = zeros(nt,nh);
para_sum = zeros(nt,nh);
integrated_total_signal = beta_atten;
integrated_single_signal = beta_atten;
int_MS = beta_atten;
beta_atten_single = beta_atten;

corrected_signal = beta_atten;
for it= 1:nt
%     lower_bound(it) = id_cb_lidar(it)+1;
%     upper_bound(it) = I_max(it)+up+5;

    perp_sum (it,lower_bound(it):upper_bound(it)) = cumsum( perp(it, ...
        lower_bound(it):upper_bound(it)),2);
    para_sum (it,lower_bound(it):upper_bound(it)) = cumsum( para(it, ...
        lower_bound(it):upper_bound(it)),2);
    
    for jz = 1:length(lower_bound(it):upper_bound(it))
    perp_int(it,lower_bound(it)+jz) = perp_sum (it,lower_bound(it)+jz)* height_res;
   
    para_int(it,lower_bound(it)+jz) = para_sum (it,lower_bound(it)+jz)* height_res;
    
    int_MS(it,lower_bound(it)+jz) = (perp_sum(it,lower_bound(it)+jz) + ...
        para_sum(it,lower_bound(it)+jz)) * height_res ;
    end 
%     depol_ratio(it,:) = perp_int(it,:) ./ para_int(it,:) ;
%     correction_factor(it,:) = (1 - depol_ratio(it,:)) ./ (1 + depol_ratio(it,:));
%     correction_factor(isnan(correction_factor)) = 1;
%     
%     int_SS = int_MS .* correction_factor ;
        
end
depol_ratio = perp_int ./ para_int ;
correction_factor = (1 - depol_ratio) ./ (1 + depol_ratio);
correction_factor(isnan(correction_factor)) = 1 ;
int_SS = int_MS .* correction_factor ;

for it= 1:nt
   for jjz = 1:length(lower_bound(it):upper_bound(it))
    beta_atten_single(it,lower_bound(it)+jjz) = (int_SS(it,lower_bound(it)+jjz) - ...
        int_SS(it,lower_bound(it)+jjz-1)) ./ height_res ; 
   end 
end
% % 
% % % corrected
% % 
% % corrected_signal = correction_factor .* beta_atten ;
sd_corrected_signal = 0.2 * corrected_signal ;

plot(beta_atten_single(50,:),height,'r',beta_atten(50,:),height,'k')

%% Aerosol Extinction-to-Backscatter ratio used in the inversion process [sr]
% and uncertainity in S used to assess the associated bias in the results 
% extinction-to-backscatter ratio for clouds should be 20
S = zeros(nt,nh);
S(:,lower_bound:upper_bound) = 20.0 ;
S(:,1:lower_bound-1) = 50.0 ;
% S = 20.0 ;     
dS = 20.0 ; 
 
% %% Prelocating matrices
% %% Prelocating matrices
%  
%   alpha=zeros(nt,nh)           ; %Extinction (alpha)
%   alpha_u=zeros(nt,nh)         ; %Extinction (alpha) dervied using upper estimate of S
%   alpha_l=zeros(nt,nh)         ; %Extinction (alpha) dervied using lower estimate of S
%   
%   alphap=zeros(nt,nh)           ; %Derivative of the Extinction (alpha)
%   alphap_u=zeros(nt,nh)         ; %Derivative of the Extinction (alpha) dervied using upper estimate of S
%   alphap_l=zeros(nt,nh)         ; %Derivative of the Extinction (alpha) dervied using lower estimate of 
% 
%   beta=zeros(nt,nh)            ; %Backscatter 
%   beta_u=zeros(nt,nh)          ; %Backscatter dervied using upper estimate of S
%   beta_l=zeros(nt,nh)          ; %Backscatter dervied using lower estimate of S
% 
%   R=zeros(nt,nh)               ; %Scattering ratio (See Eq 22)
%   R_u=zeros(nt,nh)             ; %Scattering ratio (See Eq 22) dervied using upper estimate of S
%   R_l=zeros(nt,nh)             ; %Scattering ratio (See Eq 22) dervied using lower estimate of S
% 
%   sd_alpha=zeros(nt,nh)        ; %Noise error in derivied extinction
%   sd_beta=zeros(nt,nh)         ; %Noise error in derivied backscatter
%   sd_R=zeros(nt,nh)            ; %Noise error in derivied Scattering ratio
% 
%   sd_alpha_u=zeros(nt,nh)      ; %Noise error in derivied extinction dervied using upper estimate of S
%   sd_beta_u=zeros(nt,nh)       ; %Noise error in derivied backscatter dervied using upper estimate of S
%   sd_R_u=zeros(nt,nh)          ; %Noise error in derivied Scattering ratio dervied using upper estimate of S
% 
%   sd_alpha_l=zeros(nt,nh)      ; %Noise error in derivied extinction dervied using lower estimate of S
%   sd_beta_l=zeros(nt,nh)       ; %Noise error in derivied backscatter dervied using lower estimate of S
%   sd_R_l=zeros(nt,nh)          ; %Noise error in derivied Scattering ratio dervied using lower estimate of S
% 
%   d_alpha_bias=zeros(nt,nh)    ; %Bias error in extinction due to S uncertainty
%   d_beta_bias=zeros(nt,nh)     ; %Bias error in backscatter due to S uncertainty
%   d_R_bias=zeros(nt,nh)        ; %Bias error in backscatter uncertanity due to S uncertainty
% 
%   depol_a=zeros(nt,nh)         ; %Aerosol depolarization
% 
%   Pp=zeros(nh)                 ; %P' See Eq.(11)
%   Pp_u=zeros(nh)               ; %P' corresponding to S+dS 
%   Pp_l=zeros(nh)               ; %P' corresponding to S-dS
% 
%   sd_Pp=zeros(nh)              ; %Estimated sd of P'
%   sd_Pp_u=zeros(nh)            ; %..
%   sd_Pp_l=zeros(nh)            ; %..
%   
%   height_norm_1 = zeros(1,nt);     % peak value of the signal in each time step - lower limit
%   height_norm_2 = zeros(1,nt);     % peak value of the signal in each time step - upper limit
%   
%   alpha_o_p = zeros(1,nt) ;        % Value of normalization alpha
%   alpha_o_p_u = zeros(1,nt) ;      % Value of normalization alpha for the upper limit
%   alpha_o_p_l = zeros(1,nt) ;      % % Value of normalization alpha for the lower limit
%   
%   int_S_beta = zeros(nt,nh) ;

%% Starting program main loop in time

for it = 1:nt,
      
 %% Normalization height
 % 1. Find the peak value of the signal in each time step

 height_norm_1(it) = height(I_max(it)+down) ;
 height_norm_2(it) = height(I_max(it)+up) ;  
 
 int_S_beta = 0.0 ;
 int_S_beta_u = 0.0 ;
 int_S_beta_l = 0.0 ;

%% Calculating Pp, Pp_u, Pp_l and associated standard deviations in a loop
   for iz = 1:nh
       
    % Calculations for Pp - create a function from this       
        int_S_beta = int_S_beta + (beta_m(it,iz) * S(it,iz)) * height_res ;
        Pp_corr(iz) = S(it,iz) * beta_atten_single(it,iz)  * exp(-2.0 * int_S_beta) * exp(2 * tau_mol(it,iz));
        Pp(iz) = S(it,iz) * beta_atten(it,iz)  * exp(-2.0 * int_S_beta) * exp(2 * tau_mol(it,iz));
%         Pp(it,iz) = S * beta_atten(it,iz)  * exp(-2.0 * int_S_beta) * exp(2 * tau_mol(it,iz)); 
        sd_Pp(iz) = S(it,iz) * sd_beta_atten(it,iz) * exp(-2.0 * int_S_beta) * exp(2 * tau_mol(it,iz));
        sd_Pp_corr(iz) = S(it,iz) * sd_corrected_signal(it,iz) * exp(-2.0 * int_S_beta) * exp(2 * tau_mol(it,iz));
%      % Calculations for Pp_u (upper limit)
%         S_u = S(it,iz) + dS;
%         int_S_beta_u = int_S_beta_u + (beta_m(it,iz) * S_u) * height_res ;
%         Pp_u(it,iz) = S_u * beta_atten(it,iz)  * exp(-2.0 * int_S_beta_u) * exp(2 * tau_mol(it,iz)); 
%         sd_Pp_u(it,iz) = S_u * sd_beta_atten(it,iz) * exp(-2.0 * int_S_beta_u) * exp(2 * tau_mol(it,iz));
%       
%      % Calculations for Pp_l (lower limit)
%         S_l = S(it,iz) - dS; 
%         int_S_beta_l = int_S_beta_l + (beta_m(it,iz) * S_l) * height_res ;
%         Pp_l(it,iz) = S_l * beta_atten(it,iz)  * exp(-2.0 * int_S_beta_l) * exp(2 * tau_mol(it,iz)); 
%         sd_Pp_l(it,iz) = S_l * sd_beta_atten(it,iz) * exp(-2.0 * int_S_beta_l) * exp(2 * tau_mol(it,iz));

    end
 %% Inverting lidar signal to estimate alpha_a and beta_a 
  % using Pp, Pp_u and Pp_l
  % 1. Calculating alpha_o_p, alpha_o_p_l and alpha_o_p_u
    % Prelocate matrices
%     y = zeros(it,I_max(it)+down:I_max(it)+up);
%     poly_coeff = zeros(it,I_max(it)+down:I_max(it)+up);
%     struct_S = zeros(it,I_max(it)+down:I_max(it)+up);
%     mu = zeros(it,I_max(it)+down:I_max(it)+up);
  
  x = height(I_max(it)+down:I_max(it)+up);
  y_corr(it,:) = corrected_signal(it,I_max(it)+down:I_max(it)+up);
  y(it,:) = beta_atten(it,I_max(it)+down:I_max(it)+up);
  
 [poly_coeff(it,:) , struct_S(it,:), mu(it,:)] = polyfit(x,y(it,:),1);
 [poly_coeff_corr(it,:) , struct_S_corr(it,:), mu_corr(it,:)] = polyfit(x,y_corr(it,:),1);
 
  slope(it) = poly_coeff(it,1);
  slope_corr(it) = poly_coeff_corr(it,1);

  [yfit, d_alpha_o] = polyval(poly_coeff(it,:),x,struct_S(it,:));
  yresid = y(it,:) - yfit;
  SSresid = sum(yresid.^2);
  SStotal = (length(y)-1) * var(y);
  rsq = 1 - SSresid./SStotal;
  
  % Calculate the index for io  Find the index corresponding to zo
height_o(it) = (height_norm_1(it)+ height_norm_2(it))/2;
% zo = (zn1+zn2)/2 ;
io(it)  = find(height <= height_o(it), 1, 'last' );
  
% Estimating value of alpha_o for 3 separate inversions
  alpha_o_p(it) = -0.5 * (slope(it) * 1/mean(Pp(I_max(it)+down:I_max(it)+up))) ...
                  + mean(beta_m(it,I_max(it)+down:I_max(it)+up)) * S(it,iz);
  alpha_o_p_corr(it) = -0.5 * (slope_corr(it) * 1/mean(Pp_corr(I_max(it)+down:I_max(it)+up))) ...
                  + mean(beta_m(it,I_max(it)+down:I_max(it)+up)) * S(it,iz);
%   alpha_o_p_u(it) = -0.5 * (slope(it) * 1/mean(Pp(it,I_max(it)+down:I_max(it)+up))) ...
%                   + mean(beta_m(it,I_max(it)+down:I_max(it)+up)) * S_u;
%   alpha_o_p_l(it) = -0.5 * (slope(it) * 1/mean(Pp(it,I_max(it)+down:I_max(it)+up))) ...
%                   + mean(beta_m(it,I_max(it)+down:I_max(it)+up)) * S_l;            
%  %% 2. Call the Klett function 3 times - for S, S+dS and S-dS
%   [alphap, sd_alphap] = Klett_inversion(height(I_max(it)+down:I_max(it)+up), beta_atten(it,I_max(it)+down:I_max(it)+up), ...
%     sd_Pp(it,I_max(it)+down:I_max(it)+up), height_norm_1(it), height_norm_2(it), alpha_o_p(it)) ;
  [alphap(it,:), sd_alphap(it,:)] = Klett_inversion (height, Pp, sd_Pp, ...
      io(it), alpha_o_p(it)) ;
  [alphap_corr(it,:), sd_alphap_corr(it,:)] = Klett_inversion(height, Pp_corr, sd_Pp, ...
      io(it), alpha_o_p_corr(it)) ;
%   [alphap_u(it,:), sd_alphap_u] = Klett_inversion(height, Pp_u(it,:), sd_Pp_u(it,:), ...
%       height_norm_1(it), height_norm_2(it), alpha_o_p_u(it)) ;
%   [alphap_l(it,:), sd_alphap_l] = Klett_inversion(height, Pp_l(it,:), sd_Pp_l(it,:), ...
%       height_norm_1(it), height_norm_2(it), alpha_o_p_l(it)) ;
%   
 %% Transforming alpha_p and sd_alpha_p to alpha and sd_alpha
  alpha(it,:) = alphap(it,:) ;%- S * beta_m(it,:) ; % becuse there is no aerosols below the cloud in ECSIM simulations
  alpha_corr(it,:) = alphap_corr(it,:) ;
  sd_alpha(it,:) = sd_alphap(it,:)  ;%- S * beta_m(it,:) ;
%   alpha_u(it,:) = alphap_u(it,:) ;%- S_u * beta_m(it,:) ;
%   alpha_l(it,:) = alphap_l(it,:) ;%- S_l * beta_m(it,:) ;
   
  % Estimating bias error in alpha due to the uncertainty in S
%   d_alpha_bias(it,:) = abs((alpha_u(it,:) - alpha_l(it,:))/2) ;
   
  % Calculating the corresponding beta (backscattering) values
  beta(it,:) = alpha(it,:) ./ S(it,iz) ;
  sd_beta(it,:) = sd_alpha(it,:) ./ S(it,iz) ;
%   beta_u(it,:) = alpha_u(it,:) ./ S_u ;
%   sd_beta_u(it,:) = sd_alpha_u(it,:) ./ S_u ;
%   beta_l(it,:) = alpha_l(it,:) ./ S_l ;
%   sd_beta_l(it,:) = sd_alpha_l(it,:) ./ S_l ;
%   d_beta_bias(it,:) = abs((beta_u(it,:) - beta_l(it,:))./2) ;  
%  
  
end
       

tau = height_res * cumsum(alpha,1);         %Optical depth
sd_tau = height_res * cumsum(sd_alpha,1) ;  %Estimated random error in optical depth
% tau_u = height_res * cumsum(alpha_u,1);     %Optical depth corresponding to upper S estimate
% tau_l = height_res * cumsum(alpha_l,1);     %Optical depth corresponding to lower S estimate

figure
plot(alpha(1,:),height,'.-k')
hold on
plot(Ext(1,:), height, '.-b')
hold on 
plot(alpha_corr(1,:), height, '*-r')
%refline(0,height(1,ice),'r')
% cloud_base = refline(0,height(1,icb));
% set(cloud_base,'Color','red')
% title('Retrieved Extinction coefficient');
 