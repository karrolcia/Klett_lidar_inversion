% Lidar inversion - module for simulated data
% clear all 
format long
clear all;
% load('les_1-merged_with_driz.mat')
% load('true_les_1-merged_with_driz.mat')
load('fire_27_merged_with_mask_v7.mat')
load('fire_27_merged_with_mask_v7_ss.mat')
load('true_fire_27_merged_with_mask_v7.mat')





%% Use modelled Ext to create beta_atten and make Klett inversion on it
% Calculate rho & rayleigh scattering
% Constants used in the UV_lidar_inversion

 nh = nz-1;   
 nt = nx;
 time = linspace(1,nt,nt);
 height = height(1,1:nh) .* 1.0e3;
 height_res = height(2) - height(1);  %Height resolution in km
 
 beta_atten_ss = lid_total_ss;% .* (1./(2.*Ext.*height_res).*(1-exp(-2.*Ext.*height_res)));

 av_num = 6.0222e+23;  % molecules/mole
 rg = 8.314;  % J/mole*K universal gass cnst
 rho_atm = (P(1:nt,1:nh) .* 1.0e2) ./ (rg .* T(1:nt,1:nh)) .* av_num; % pressure from mb to Pa
 rho_atm = rho_atm .* 1.0e-6  ; %in 1/cm3

ray_sigma_beta
% Calculating molecular extinction and backscatter fields
alpha_m = rho_atm * alpha_ray ./1.0e3 ;  % now in 1/m - those values are correct (should be 1.0e-10)
beta_m = rho_atm * beta_ray   ./1.0e3 ; % now in 1/m

tau_mmol = zeros(nt,nh);
tau_cld = zeros(nt,nh);

for it = 1:nt
%     tau_mol(it,:)= cumsum(alpha_m(it,:));
%     tau_mol(it,:)=tau_mol(it,:) .* height_res;

%     tau_cld(it,:) =cumsum(Ext(it,1:nh));
%     tau_cld(it,:) = tau_cld(it,:) .* height_res;
%     int_sig = zeros(1,n);
tau_mmol(1) = 0.5 * alpha_m(it,1)* height(1)^2;
tau_cld(1) = 0.5 * Ext(it,1) * height(1)^2;

for i=2:nh-1 
   tau_mmol(it,i) = tau_mmol(it,i-1) + 0.5*((alpha_m(it,i-1)+alpha_m(it,i))*height_res); 
   tau_cld(it,i) = tau_cld(it,i-1) +  0.5*((Ext(it,i-1)+ Ext(it,i))*height_res);
   
end
end

alpha_total = Ext ;


S = zeros (nt,nh);
S(:, :) = 20;
% S(:,1:13) = 20;
beta_atten = ((alpha_total./S)+beta_m) .* exp(-2.0 .* (tau_cld+tau_mmol)) ;
beta_atten_org = beta_atten;


sd_beta_atten = 0.2 .* beta_atten ;

 % Conver Range Corrected Signal to Power
P_sig = zeros(nt,nh);
P_sig_org = zeros(nt,nh);
P_sig_noise = zeros(nt,nh);
noise = zeros(nt,nh);
c = 3.0e8;
planck = 6.62607004 * 1.0e-34; % Planck's Constant
wave = 355 * 1.0e-9 ; % lidar wavelength
mirror_dia = 0.45; % Receiver mirror radius
pulse_power = 0.35; % based on UV Raman from CHilbolton (in J)
pulse_length = 7 *1.0e-9 ;% based on UV Raman from CHilbolton (in s) - original 7 ns
dR = height_res; 
efficiency = 0.8;
N_l = (pulse_power*pulse_length)/((planck*c)/wave)*3000; % times 3000 because otherwise it's calculated for a single pulse
C_lid =  N_l*dR*((mirror_dia/2)^2*pi)*efficiency;
N_b = 20;
for it = 2:nt
for jz = 1:nh
 P_sig_org(it,jz) = beta_atten(it,jz) ./ height(jz).^2;
% %  C_lid(it,:) = (height.^2 .* P_sig_org(it,:)) ./ rho_atm(it,:);
photons(it,jz) = (P_sig_org(it,jz)).*(C_lid);
noise(it,jz) = (sqrt((photons(it,jz) + (N_b^2))).*randn())/C_lid;%.*abs(randn());
P_sig_noise(it,jz) = P_sig_org(it,jz) +  noise(it,jz);

snr_P_sig(it,jz) = snr(P_sig_org(it,jz),noise(it,jz));

beta_atten(it,jz) = (P_sig_noise(it,jz)) .* height(jz).^2;
end 
end
 sd_P_sig = 0.2 .* P_sig;
%  


% P_sig_noise = P_sig + sqrt((P_sig.*C_lid) + (sqrt( P_sig./C_lid).^2)).*abs(randn());
P_sig = P_sig_noise;
% figure
% plot(P_sig_org(10,:),height,'k', P_sig_noise(10,:), height,'r')


 para = beta_atten ./ (lin_depol(1:nt,1:nh) + 1) ;
 perp = beta_atten - para;
 
% Calculate LWP - convert cm to g/m^2
lwp = LWP ./ 1.e-4 ;
lwp = lwp(:,1)';

 %% Calculate normalization value for each time step
% Find the peak of the attenuated backscatter
 [Peak_P, I_max] = max(beta_atten_ss, [], 2) ;
 [Peak_perp, perp_max] = max(perp, [], 2) ;


% Aerosol Extinction-to-Backscatter ratio used in the inversion process [sr]
% and uncertainity in S used to assess the associated bias in the results 
% extinction-to-backscatter ratio for clouds should be 20
S = zeros(nt,nh);
 for it = 1:nt
  id_cb_lidar(it) = find (perp(it,1:perp_max(it)) >= (Peak_perp(it)/10), 1,'first'); 
  id_cb_radar(it) = find((Ze(it,:)) >= 1.e-3, 1, 'first');
  id_ct_radar(it) = find((Ze(it,:)) >= 1.e-3, 1, 'last');

  norm_down(it) = I_max(it)+4;
   norm_up(it) = norm_down(it) +3 ;

  height_norm_1(it) = height(norm_down(it)) ;
  height_norm_2(it) = height(norm_up(it)) ; 
  height_o(it)= (height_norm_1(it) + height_norm_2(it))/2;
  io(it)  = find(height <= height_o(it), 1, 'last' );
  S(it, 1:id_cb_lidar(it)-1) = 50.0;
  S(it,id_cb_lidar(it):end) = 20.0 ;
%   snr_norm(it,:) = snr_P_sig(it,id_cb_lidar(it):norm_down(it) +3);

 end
dS = 20.0 ; 


 tau_mol = zeros(nt,nh);
 int_S_beta = zeros(nt,nh) ;
    
 
for it = 1:nt
      
 tau_mol(it,1) = alpha_m(it,1) * height(1)^2;
 int_S_beta(it,1) = (beta_m(it,1) * S(it,1)) * height(1)^2;

%% Calculating Pp, Pp_u, Pp_l and associated standard deviations in a loop
 for iz = 2:nh
%       tau_cld(it,i) = tau_cld(it,i-1) +  0.5*((Ext(it,i-1)+ Ext(it,i))*height_res);
tau_mol(it, iz) = tau_mol(it,iz-1) + 0.5 * ((alpha_m(it,iz-1) + alpha_m(it,iz)) * height_res);
int_S_beta(it,iz) =  int_S_beta(it,iz-1) + 0.5 *((beta_m(it,iz-1)*S(it,iz-1) + beta_m(it,iz)*S(it,iz)) * height_res);

Pp(it,iz) = S(it,iz) * P_sig(it,iz) * exp(-2.0*int_S_beta(it,iz)) * exp(2*tau_mol(it,iz));
sd_Pp(it,iz) = S(it,iz) * sd_P_sig(it,iz) * exp(-2.0 *int_S_beta(it,iz)) * exp(2*tau_mol(it,iz));

 end

 Pp_int{it} = Pp(it,:) ;
 
% Calculate the index for io  Find the index corresponding to zo
height_o(it)= (height_norm_1(it) + height_norm_2(it))/2;
io(it)  = find(height <= height_o(it), 1, 'last' );

% Calculating slope of the signal
x = log(beta_atten(it,norm_down(it):norm_up(it)));
y = height;

beta_slope{it} = x;
height_slope{it} = y;
  
alpha_o_p(it,:) = -0.5 * (diff(beta_slope{it}) ./ height_res);% +  beta_m(it,io(it)) * S(it,io(it));
alpha_slope(it) = alpha_o_p(it,2); 
% alpha_slope(it) = mean(alpha_o_p(it,:));
  
% Estimating value of alpha_o 

alpha_o(it) = Ext(it,io(it)) + S(it, io(it)) .* beta_m(it, io(it));

[alphap{it}, sd_alphap{it}] = Klett_inversion(height, Pp_int{it}, sd_Pp, ...
     io(it), alpha_slope(it)) ;

alpha(it,:) = alphap{it}  - S(it,:) .* beta_m(it,:) ; % becuse there is no aerosols below the cloud in ECSIM simulations
sd_alpha(it,:) = sd_alphap{it} ;%- S(it,:) .* beta_m(it,:) ;

% Calculating the corresponding beta (backscattering) values
beta(it,:) = alpha(it,:) ./ S(it,iz) ;
sd_beta(it,:) = sd_alpha(it,:) ./ S(it,iz) ;
  
clear Pp sd_Pp x y
end
tau_org = height_res * cumsum(Ext,1); % Modelled ECSIM optical depth       
tau = height_res * cumsum(alpha,1);         %Optical depth
sd_tau = height_res * cumsum(sd_alpha,1) ;  %Estimated random error in optical depth
e = std(alpha);

tt=20;

figure('NumberTitle','off', ...
    'Units', 'centimeters','Position',[2 20 15 15])
% errorbar(alpha(100,:),e)
plot(alpha(tt,id_cb_lidar(tt):norm_up(tt)), tau_org(tt,id_cb_lidar(tt):norm_up(tt)), '-.r')
% herrorbar(alpha(nt,:),sd_alpha(nt,:))
hold on
plot(Ext(tt,id_cb_lidar(tt):norm_up(tt)), tau_org(tt,id_cb_lidar(tt):norm_up(tt)), '-.b')
hold on 
xlabel('\alpha [m^{-1} sr]')
ylabel('\tau ')
legend('Retrieved extinction', '"True" extinction')
% normalisation = patch([0.25 -0.1 -0.1 0.25], ...
%     [height(norm_range_1) height(norm_range_1) height(norm_range_2) height(norm_range_2)], ...
%     'b');
% set(normalisation,'FaceColor','b', ...
%     'FaceAlpha', 0.1)
% % refline([0 height(I_max(tt))])
% refline([0 tau_org(norm_down)])
hold on
% plot(alpha_corr(1,:), height, '*-r')
%refline(0,height(1,ice),'r')
% cloud_base = refline(0,height(1,icb));
% set(cloud_base,'Color','red')
% title('Retrieved Extinction coefficient');

% figure
% plot(beta_atten(1,:).*200,height,'.-r' , ...
%    Ext(1,:),height,'.-b', ...
%    alpha(1,:),height,'.-g', ...
%    'LineWidth', 1.2)
% hold on
% refline([0 height(io(1))])

%% Simple check of accuracy of alpha_o
% Ext = Ext + S.*beta_m;
for i = 1:nt
    accuracy(i) = alpha_slope(i)./Ext(i,io(i));
    accuracy2(i) = (alpha_slope(i) - Ext(i,io(i)))./mean(Ext(:,io(i)));
    bias(i) = alpha_slope(i) - alpha_o(i);
end
mean_bias_alpha_o = mean(bias)
mean_accuracy_alpha_o = mean(accuracy)
mean_accuracy_2 = mean(accuracy2)

%% Calculate statistical data and make tables

for it = 1:nt
retrieved_cloud{it} = id_cb_lidar(it):norm_down(it);
height_cloud{it} = height(retrieved_cloud{it});
end


for it = 1:nt
    for iih = 1:nh
      accu_alpha_plot(it,iih) = (alpha(it,iih) ./ Ext(it,iih)); 
    end 
    
    for ih = 1:length(retrieved_cloud{it})
accu_alpha(it,ih) = (alpha(it,retrieved_cloud{it}(ih)) ./ Ext(it,retrieved_cloud{it}(ih)));
err_alpha(it,ih) = (alpha(it,retrieved_cloud{it}(ih)) - Ext(it,retrieved_cloud{it}(ih))); 
RMSE(:,ih) = sqrt(mean((Ext(it,retrieved_cloud{it}(ih)) - alpha(it,retrieved_cloud{it}(ih))).^2));
bias_alpha(:,ih) = (nansum(alpha(:,retrieved_cloud{it}(ih))) - nansum(Ext(:,retrieved_cloud{it}(ih)))) ./ nt;

    end  
accu_alpha_plot(it,1:retrieved_cloud{it}-1)=NaN;
accu_alpha_plot(it,max(retrieved_cloud{it})+1:end)=NaN;
    
end
accu_alpha_plot(accu_alpha_plot ==0) = NaN;
accu_alpha(accu_alpha ==0) = NaN;
err_alpha(err_alpha==0)= NaN;
% bias_alpha(bias_alpha==0)= NaN;
%% Latex table
fprintf('height above cb & mean accu_alpha & mean alpha_error & bias & RMSE \n')

for ii=1:6
    fprintf('%8.0f & %8.12f &  %8.12f &  %8.12f &  %8.12f &   \\\\ \n', ...
        height(retrieved_cloud{ii}(ii)) - height(retrieved_cloud{ii}(1)) , ...
        nanmean(accu_alpha(:,ii)), ...
        nanmean(err_alpha(:,ii)), ...
        bias_alpha(ii), ...
        mean(RMSE(ii)))
end

figure
plot(alpha(15,1:50),height(1:50),'.-k',Ext(15,1:50),height(1:50),'.-r')
hold on
refline([0 height(id_cb_lidar(1))])
refline([0 height(norm_down(1))])
refline([0 height(io(1))])


% beta_atten_org2 = beta_atten_org .* (1./(2.*alpha.*height_res).*(1-exp(-2.*alpha.*height_res)));

figure
semilogx(log(beta_atten_ss(100,1:50)),height(1:50),'.-r',log(beta_atten_org(100,1:50)),height(1:50),'.-b')
legend(['ECSIM ATB'], ['Calculated ATB'])
% log(beta_atten_ss(1,:))-log(beta_atten_org(1,:))