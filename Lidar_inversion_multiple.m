% Lidar inversion - module for simulated data
clear all 
format long
interp_factor = 5; % vertical resolution in the interpolated signal
% % 
load('les_1-merged_with_driz_5m.mat')
load('true_les_1-merged_with_driz_5m.mat')
% load('fire_27_merged_with_mask_v6.mat')
% load('fire_27_merged_with_mask_v6_ss.mat')
% load('true_fire_27_merged_with_mask_v6.mat')


beta_atten = lid_total ; % Original unit is 1/m/sr, change to 1/km/sr

sd_beta_atten = 0.2 .* beta_atten ;

 nh = nz;   
 nt = nx;
 time = linspace(1,nt,nt);
 height = height(1,:) .* 1.0e3;
 height_res = height(2) - height(1);  %Height resolution in km
 
%  height_int = min(height):5:max(height);
  height_int = interp(height,interp_factor);

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
 alpha_m = rho_atm .* alpha_ray ./1.0e3 ;  % in 1/km chnaged to 1/m
 beta_m = rho_atm .* beta_ray  ./1.0e3 ; % 
for it = 1:nt 
alpha_m_int(it,:) = interp(alpha_m(it,:),interp_factor);
beta_m_int(it,:) = interp(beta_m(it,:),interp_factor);
end  
%% Find peak signal value in each time step

% Find the peak of the attenuated backscatter
 [Peak_P, I_max] = max(beta_atten, [], 2) ;
 [Peak_perp, perp_max] = max(perp, [], 2) ;
 
 %% Calculate normalization value for each time step
% Define up and down limit for normalization
down = 4;
up = 7;

 for it = 1:nt
    id_cb_lidar(it) = find (perp(it,1:perp_max(it)) >= (Peak_perp(it)/10), 1,'first');   
    id_cb_radar(it) = find((Ze(it,:)) >= 1.e-3, 1, 'first');
    id_ct_radar(it) = find((Ze(it,:)) >= 1.e-3, 1, 'last');
    norm_down(it) = find(beta_atten(it,:) <=  1.e-6, 1, 'first');
    norm_up(it) = norm_down(it) +3 ;
% Calculate the index for io  Find the index corresponding to zo
%     norm_range_1 = I_max(it)  + down;
%     norm_range_2 = I_max(it)  + up;
    height_norm_1(it) = height(norm_down(it));
    height_norm_2(it) = height(norm_up(it)) ;
%     height_norm_1(it) = height(norm_range_1) ;
%     height_norm_2(it) = height(norm_range_2) ; 
%     height_o(it)= (height_norm_1(it) + height_norm_2(it))/2;
    height_o(it) = (height_norm_1(it) + height_norm_2(it))/2;
    io(it)  = find(height <= height_o(it), 1, 'last' );
    io_int(it)  = find(height_int <= height_o(it), 1, 'last' );
 end

%% Corecting signal for the multiple scattering
perp_int = zeros(nt,nh);
para_int = zeros(nt,nh);

lower_bound = id_cb_lidar-1;
upper_bound = id_ct_radar;
perp_sum = zeros(nt,nh);
para_sum = zeros(nt,nh);
% integrated_total_signal = beta_atten;
% integrated_single_signal = beta_atten;
int_MS = beta_atten;
beta_atten_single = beta_atten;
P_sig_single = P_sig ;

corrected_signal = beta_atten;
for it= 1:nt

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
        
end
depol_ratio = perp_int ./ para_int ;
correction_factor = (1 - depol_ratio) ./ (1 + depol_ratio);
correction_factor(isnan(correction_factor)) = 1 ;
int_SS = int_MS .* correction_factor ;

for it= 1:nt
   for jjz = 1:length(lower_bound(it):upper_bound(it))
    beta_atten_single(it,lower_bound(it)+jjz) = (int_SS(it,lower_bound(it)+jjz) - ...
        int_SS(it,lower_bound(it)+jjz-1)) ./ height_res ; 
    P_sig_single(it,lower_bound(it)+jjz) = beta_atten_single(it,lower_bound(it)+jjz) ...
        ./ (height(lower_bound(it)+jjz).^2) ;
   end 
end

sd_P_sig_single = 0.2 * corrected_signal ;

 P_log = log(P_sig_single);
%  P_sig_int = interp(P_log,interp_factor);
%  sd_P_sig_int = interp(sd_P_sig_single,interp_factor);
 for ii = 1:nt
 P_sig_int(ii,:) = interp(P_log(ii,:),interp_factor);    
%  P_sig_int(ii,:) = interp1(height,P_sig_single(ii,:),height_int,'linear');
 sd_P_sig_int(ii,:) = interp1(height,sd_P_sig_single(ii,:),height_int,'linear'); 
 end
 P_sig_int = exp(P_sig_int);

% figure('NumberTitle','off', ...
%     'Units', 'centimeters','Position',[2 20 15 15])
% plot(beta_atten_single(50,:),height,'r',beta_atten(50,:),height,'k')
% 
% figure('NumberTitle','off', ...
%     'Units', 'centimeters','Position',[2 20 15 15])
% plot(P_sig(50,:),height,'r',P_sig_single(50,:),height,'k',P_sig_int(50,:),height_int,'b')

%% Aerosol Extinction-to-Backscatter ratio used in the inversion process [sr]
% and uncertainity in S used to assess the associated bias in the results 
% extinction-to-backscatter ratio for clouds should be 20
S = zeros(nt,nh);
S(:,lower_bound:end) = 20.0 ;
S(:,1:lower_bound-1) = 50.0 ;
% S(:,:) = 20.0;
S_int =  zeros(nt,length(height_int));
S_int(:,lower_bound:end) = 20.0 ;
S_int(:,1:lower_bound-1) = 50.0 ;
% S_int(:,:) = 20.0;
% % S = 20.0 ;     
dS = 20.0 ; 
 
for it = 1:nt
      

 tau_mol = 0;
 int_S_beta = 0 ; 

%% Calculating Pp, Pp_u, Pp_l and associated standard deviations in a loop
 for iz = 1:nh
      
tau_mol = tau_mol + alpha_m(it,iz) * height_res;
int_S_beta = int_S_beta + (beta_m(it,iz) * S(it,iz)) * height_res ;

Pp(iz) = S(it,iz) * P_sig(it,iz) * exp(-2.0*int_S_beta) * exp(2*tau_mol);
Pp_single(iz) = S(it,iz) * P_sig_single(it,iz) * exp(-2.0*int_S_beta) * exp(2*tau_mol);
  
sd_Pp(iz) = S(it,iz) * sd_P_sig(it,iz) * exp(-2.0 *int_S_beta) * exp(2*tau_mol);
sd_Pp_single(iz) = S(it,iz) * sd_P_sig_single(it,iz) * exp(-2.0 *int_S_beta) * exp(2*tau_mol);

 end

 Pp_int{it} = Pp ;
 Pp_int_s{it} = Pp_single;

 

 for izn = 1:length(height_int)
  tau_mol_int = 0;
 int_S_beta_n = 0 ;      
tau_mol_int = tau_mol_int + alpha_m_int(it,izn) * height_res;
int_S_beta_n = int_S_beta_n + (beta_m_int(it,izn) * S_int(it,izn)) * height_res ;

Pp_n(izn) = S_int(it,izn) * P_sig_int(it,izn) * exp(-2.0*int_S_beta_n) * exp(2*tau_mol_int);
sd_Pp_interp(izn) = S_int(it,izn) * sd_P_sig_int(it,izn) * exp(-2.0 *int_S_beta_n) * exp(2*tau_mol_int);

 end

 Pp_interp{it} = Pp_n ;
 
 
% Calculate the index for io  Find the index corresponding to zo
height_o(it)= (height_norm_1(it) + height_norm_2(it))/2;
io(it)  = find(height <= height_o(it), 1, 'last' );

io_int(it)  = find(height_int <= height_o(it), 1, 'last' );

% Calculating slope of the signal
% Calculating slope of the signal  
% x = log(beta_atten(it,I_max(it)+down:I_max(it)+up));
% y = height(I_max(it)+up) - height(I_max(it)+down);
% x_s = log(beta_atten_single(it,I_max(it)+down:I_max(it)+up));

% x = log(beta_atten(it,(norm_range_1:norm_range_2)));
% y = height(norm_range_1:norm_range_2);
% x_s = log(beta_atten_single(it,(norm_range_1:norm_range_2)));

x = log(beta_atten(it,(norm_down:norm_up)));
y = height(norm_down:norm_up);
x_s = log(beta_atten_single(it,(norm_down:norm_up)));

beta_slope{it} = x;
beta_slope_s{it} = x;
height_slope{it} = y;
  
% [poly_coeff , struct_S, mu] = polyfit(beta_slope{it},height_slope{it},1);
% slope = poly_coeff(1);

%   [yfit, d_alpha_o] = polyval(poly_coeff,x,struct_S);
%   yresid = y - yfit;
%   SSresid = sum(yresid.^2);
%   SStotal = (length(y)-1) * var(y);
%   rsq = 1 - SSresid./SStotal;

  
%   lm = fitlm(x,y);
alpha_o_p(it,:) = -0.5 * (diff(beta_slope{it}) ./ height_res);% +  beta_m(it,io(it)) * S(it,io(it));
% alpha_slope(it) = alpha_o_p(it,round(length(beta_slope{it})/2)); 
alpha_slope(it) = mean(alpha_o_p(it,:))*1.2;

alpha_o_p_s(it,:) = -0.5 * (diff(beta_slope_s{it}) ./ height_res);% +  beta_m(it,io(it)) * S(it,io(it));
% alpha_slope_s(it) = alpha_o_p_s(it,round(length(beta_slope_s{it})/2)); 
alpha_slope_s(it) = mean(alpha_o_p_s(it,:))*1.2;

% alpha_slope(it) = alpha_o_p(it,1); 
% The io is in the middle of the normalisation interval
% Calculate RCS_o for Klett inversion - mean values of range corrected
% signal over the normalisation interval
% RCS_o(it) = mean((Pp_int{it}(I_max(it)+down:I_max(it)+up)).*height(I_max(it)+down:I_max(it)+up).^2)
  
% Estimating value of alpha_o 

alpha_o(it) = Ext(it,io(it)) + S(it, io(it)) .* beta_m(it, io(it));
% alpha_o_p(it) = 0.2;
%  alpha_o_p(it) = -0.5 * (slope(it));% + (1/(Pp_int{it}(io(it)))+ beta_m(it,io(it)) * S(it,iz)));

[alphap{it}, sd_alphap{it}] = Klett_inversion(height, Pp_int{it}, sd_Pp, ...
     io(it), alpha_slope(it)) ;
[alphap_s{it}, sd_alphap_s{it}] = Klett_inversion(height, Pp_int_s{it}, sd_Pp_single, ...
     io(it), alpha_slope(it)) ; 
 
[alphap_interp{it}, sd_alphap{it}] = Klett_inversion(height_int, Pp_interp{it}, sd_Pp_interp, ...
     io_int(it), alpha_slope_s(it)) ;
 


alpha(it,:) = alphap{it} - S(it,:) .* beta_m(it,:) ; % becuse there is no aerosols below the cloud in ECSIM simulations
sd_alpha(it,:) = sd_alphap{it} ;%- S(it,:) .* beta_m(it,:) ;

alpha_s(it,:) = alphap_s{it}- S(it,:) .* beta_m(it,:)   ;%-
alpha_int(it,:) = alphap_interp{it};

   
 % Calculating the corresponding beta (backscattering) values
beta(it,:) = alpha(it,:) ./ S(it,iz) ;
sd_beta(it,:) = sd_alpha(it,:) ./ S(it,iz) ;
  
clear Pp sd_Pp x y
end
       
tau_org = height_res * cumsum(Ext,1);
tau = height_res * cumsum(alpha,1);         %Optical depth
sd_tau = height_res * cumsum(sd_alpha,1) ;  %Estimated random error in optical depth
e = std(alpha);

tt=7;

figure('NumberTitle','off', ...
    'Units', 'centimeters','Position',[2 20 15 15])
plot(alpha(tt,:), tau_org(tt,:), '.r','MarkerSize',12)
% herrorbar(alpha(nt,:),sd_alpha(nt,:))
hold on
plot(Ext(tt,:), tau_org(tt,:), '.b','MarkerSize',12)
hold on 
plot(alpha_s(tt,:), tau_org(tt,:), '.k','MarkerSize',12)
hold on 
ylim([0 15])
% normalisation = patch([0.25 -0.1 -0.1 0.25], ...
%     [height(norm_range_1) height(norm_range_1) height(norm_range_2) height(norm_range_2)], ...
%     'b');
% set(normalisation,'FaceColor','b', ...
%     'FaceAlpha', 0.1)
% refline([0 height(id_cb_lidar(tt))])
hold on
% title('Retrieved Extinction coefficient');

%% Simple check of accuracy of alpha_o
for i = 1:nt
    accuracy(i) = alpha_slope(i)./Ext(i,io(i));
end
max(accuracy(~isinf(accuracy)))
min(accuracy(~isinf(accuracy)))
mean(accuracy(~isinf(accuracy)))

%% Calculate statistical data and make tables


%% Calculate statistical data and make tables
% retrieved_cloud = id_cb_lidar:I_max+down;
% accu_alpha = alpha(:,retrieved_cloud) ./ Ext(:,retrieved_cloud);
% accu_alpha(accu_alpha == Inf) = NaN;
% accu_alpha(accu_alpha == -Inf) = NaN;
% err_alpha = abs((alpha(:,retrieved_cloud) - Ext(:,retrieved_cloud))   ... 
%             ./ Ext(:,retrieved_cloud)) * 100;
% err_alpha(err_alpha == Inf) = NaN;      
% err_alpha(err_alpha == -Inf) = NaN; 


for it = 1:nt
retrieved_cloud{it} = id_cb_lidar(it):I_max(it)+down;
end


for it = 1:nt
    for ih = 1:length(retrieved_cloud{it})
accu_alpha(it,ih) = (alpha(it,retrieved_cloud{it}(ih)) ./ Ext(it,retrieved_cloud{it}(ih)));
err_alpha(it,ih) = (abs((alpha(it,retrieved_cloud{it}(ih)) - Ext(it,retrieved_cloud{it}(ih)))   ... 
            ./ Ext(it,retrieved_cloud{it}(ih))) * 100); 
accu_alpha_sscat(it,ih) = (alpha_s(it,retrieved_cloud{it}(ih)) ./ Ext(it,retrieved_cloud{it}(ih)));
err_alpha_sscat(it,ih) = (abs((alpha_s(it,retrieved_cloud{it}(ih)) - Ext(it,retrieved_cloud{it}(ih)))   ... 
            ./ Ext(it,retrieved_cloud{it}(ih))) * 100);  
 single_improvement(it,ih) =  ((err_alpha(it,ih) - err_alpha_sscat(it,ih))) ;   
 bias_alpha(:,ih) = (sum(alpha(:,retrieved_cloud{it}(ih))) - sum(Ext(:,retrieved_cloud{it}(ih)))) ./ nt;
    end
   
end

%% Latex table
fprintf('height above cb & mean alpha_error & max alpha_error & min alpha_error & mean alpha_error_single & max alpha_error_single & min alpha_error & single improvement \n')
for ii=1:length(retrieved_cloud{1})
    fprintf('%8.0f & %8.2f &  %8.2f &  %8.2f & %8.2f &  %8.2f &  %8.2f %8.2f \\\\ \n', ...
        height(retrieved_cloud{ii}(ii)) - height(retrieved_cloud{ii}(1)) , ...
        mean(err_alpha(:,ii)), ...
        max(err_alpha(:,ii)), ...
        min(err_alpha(:,ii)),...
        mean(err_alpha_sscat(:,ii)), ...
        max(err_alpha_sscat(:,ii)), ...
        min(err_alpha_sscat(:,ii)), ...
        mean(single_improvement(:,ii)))

end
 