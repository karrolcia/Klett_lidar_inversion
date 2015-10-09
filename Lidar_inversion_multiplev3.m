% Lidar inversion - module for simulated data
clear all 
format long

load('fire_27_merged_with_mask_v7.mat')
load('fire_27_merged_with_mask_v7_ss.mat')
load('true_fire_27_merged_with_mask_v7.mat')

bias = 0;

beta_atten = lid_total; % Original unit is 1/m/sr, change to 1/km/sr
% s = RandStream('mt19937ar','seed',5489);
% beta2 = awgn(beta_atten, 5, 0, s);
% beta_atten = beta2;
sd_beta_atten = 0.2 .* beta_atten ;

 nh = nz-1;   
 nt = nx;
 time = linspace(1,nt,nt);
 height = height(1,:) .* 1.0e3;
 height_res = height(2) - height(1);  %Height resolution in km


% Conver Range Corrected Signal to Power & addnoise 

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
N_b = 200;
for it = 1:nt
for jz = 1:nh

P_sig_org(it,jz) = beta_atten(it,jz) ./ height(jz).^2;
photons(it,jz) = (P_sig_org(it,jz)).*(C_lid);
noise(it,jz) = (sqrt((photons(it,jz) + (N_b^2))).*randn())/C_lid;%.*abs(randn());
P_sig_noise(it,jz) = P_sig_org(it,jz) +  noise(it,jz);

snr_P_sig(it,jz) = snr(P_sig_org(it,jz),noise(it,jz));

beta_atten_n(it,jz) = (P_sig_noise(it,jz)) .* height(jz).^2;
end 
end

P_sig = P_sig_noise;
sd_P_sig = 0.2 .* P_sig;
 
 para = beta_atten_n ./ (lin_depol + 1) ;
 perp = beta_atten_n - para;
 
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

 %% Find peak signal value in each time step

% Find the peak of the attenuated backscatter
 [Peak_P, I_max] = max(beta_atten, [], 2) ;
 [Peak_perp, perp_max] = max(perp, [], 2) ;
 
 for it = 1:nt
    id_cb_lidar(it) = find (perp(it,1:perp_max(it)) >= (Peak_perp(it)/10), 1,'first'); 
%     id_cb_lidar(it) = find(Ext(it,:) > 0,1, 'first');
    id_cb_radar(it) = find((Ze(it,:)) >= 1.e-3, 1, 'first');
    id_ct_radar(it) = find((Ze(it,:)) >= 1.e-3, 1, 'last');

    norm_down(it) = I_max(it)+5;
    norm_up(it) = norm_down(it) +3 ;

    height_norm_1(it) = height(norm_down(it));
    height_norm_2(it) = height(norm_up(it)) ;

    height_o(it) = (height_norm_1(it) + height_norm_2(it))/2;
    io(it)  = find(height <= height_o(it), 1, 'last' );
    
    S(it, 1:id_cb_lidar(it)-1) = 50.0;
    S(it,id_cb_lidar(it):end) = 20.0 ;

 end
%% Apply correction for the attenuation within the gridcell
%1/(2*alpha(z)*dz)*(1-exp(-2*alpha(z)*dz))
 
%% Corecting signal for the multiple scattering
perp_int = zeros(nt,nh);
para_int = zeros(nt,nh);
perp_sum = zeros(nt,nh);
para_sum = zeros(nt,nh);

int_MS = beta_atten;
beta_atten_single = beta_atten;
P_sig_single = P_sig ;
corrected_signal = beta_atten;

% for it= 1:nt
% % 
% %     perp_int (it,id_cb_lidar(it)) = 0.5 * perp(it,id_cb_lidar(it));
% %     para_int (it,id_cb_lidar(it)) = 0.5 * para(it,id_cb_lidar(it));
%     
%         perp_int (it,id_cb_lidar(it):id_ct_radar(it)) = cumsum( perp(it, ...
%         id_cb_lidar(it):id_ct_radar(it)),2).*height_res;
%     para_int (it,id_cb_lidar(it):id_ct_radar(it)) = cumsum( para(it, ...
%         id_cb_lidar(it):id_ct_radar(it)),2).*height_res;
%     int_MS(it,id_cb_lidar(it):id_ct_radar(it)) = 
%     
%     for jz = 1:length(id_cb_lidar(it):id_ct_radar(it))
%         
% %         perp_int(it,id_cb_lidar(it)+jz) = 
%         
% % %            
% %     perp_int(it,id_cb_lidar(it)+jz) = perp_int(it,id_cb_lidar(it)+jz-1) + ...
% %         0.5 * ((perp(it,id_cb_lidar(it)+jz-1)+ ...
% %         perp(it,id_cb_lidar(it)+jz))*height_res)  ;
% %    
% %     para_int(it,id_cb_lidar(it)+jz) = para_int(it,id_cb_lidar(it)+jz-1) + ...
% %         0.5 * ((para(it,id_cb_lidar(it)+jz-1) + ...
% %         para(it,id_cb_lidar(it)+jz))*height_res)  ;   
%     
%     int_MS(it,id_cb_lidar(it)+jz) = perp_int(it,id_cb_lidar(it)+jz) + ...
%         para_int(it,id_cb_lidar(it)+jz); 
%     
%     end    
% end
% depol_ratio = perp_int ./ para_int ;
% correction_factor = (1 - depol_ratio) ./ (1 + depol_ratio);
% correction_factor(isnan(correction_factor)) = 1 ;
% int_SS = int_MS .* correction_factor ;
% 
% beta_atten_single = int_SS ./ height_res;
% for it = 1:nt
% P_sig_single(it,:) = beta_atten_single(it,:) ./ height.^2;
% end
% % for it= 1:nt
% %    for jjz = 1:length(id_cb_lidar(it)-1:id_ct_radar(it))
% %     beta_atten_single(it,id_cb_lidar(it)+jjz) = (int_SS(it,id_cb_lidar(it)+jjz)) ./ height_res ; 
% %     P_sig_single(it,id_cb_lidar(it)) = beta_atten_single(it,id_cb_lidar(it)) ...
% %         ./ (height(id_cb_lidar(it)).^2) ;
% %    end 
% %    
% % end
% 
% sd_P_sig_single = 0.2 * corrected_signal ;

for it= 1:nt

    perp_sum (it,id_cb_lidar(it):id_ct_radar(it)) = cumsum( perp(it, ...
        id_cb_lidar(it):id_ct_radar(it)),2);
    para_sum (it,id_cb_lidar(it):id_ct_radar(it)) = cumsum( para(it, ...
        id_cb_lidar(it):id_ct_radar(it)),2);
    
    perp_int(it,id_cb_lidar(it))= perp_sum(it,id_cb_lidar(it))* height_res;
    para_int(it,id_cb_lidar(it)) = para_sum(it,id_cb_lidar(it))* height_res;
    
    int_MS(it,id_cb_lidar(it)) = (perp_sum(it,id_cb_lidar(it)) + ...
        para_sum(it,id_cb_lidar(it))) * height_res ;
    
    for jz = 1:length(id_cb_lidar(it):id_ct_radar(it))
    perp_int(it,id_cb_lidar(it)+jz) = perp_sum(it,id_cb_lidar(it)+jz)* height_res;
   
    para_int(it,id_cb_lidar(it)+jz) = para_sum(it,id_cb_lidar(it)+jz)* height_res;
    
    int_MS(it,id_cb_lidar(it)+jz) = (perp_sum(it,id_cb_lidar(it)+jz) + ...
        para_sum(it,id_cb_lidar(it)+jz)) * height_res ;
    end 
        
end
depol_ratio = perp_int ./ para_int ;
correction_factor = (1 - depol_ratio) ./ (1 + depol_ratio);
correction_factor(isnan(correction_factor)) = 1 ;
int_SS = int_MS .* correction_factor ;

for it= 1:nt
  beta_atten_single(it,id_cb_lidar(it)) = (int_SS(it,id_cb_lidar(it)) - ...
     int_SS(it,id_cb_lidar(it)-1)) ./ height_res;
  P_sig_single(it,id_cb_lidar(it)) = beta_atten_single(it,id_cb_lidar(it)) ./ ...
      (height(id_cb_lidar(it)).^2);
 for jjz = 1:length(id_cb_lidar(it):id_ct_radar(it))
  beta_atten_single(it,id_cb_lidar(it)+jjz) = (int_SS(it,id_cb_lidar(it)+jjz) - ...
     int_SS(it,id_cb_lidar(it)+jjz-1)) ./ height_res ;    
  P_sig_single(it,id_cb_lidar(it)+jjz) = beta_atten_single(it,id_cb_lidar(it)+jjz) ...
        ./ (height(id_cb_lidar(it)+jjz).^2) ;
 end 
   
end
sd_P_sig_single = 0.2 * corrected_signal ;

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
for it=1:nt
    S(it,id_cb_lidar(it):end) = 20.0 ;
    S(it,1:id_cb_lidar(it)-1) = 50.0 ;
end
dS = 20.0 ; 

 tau_mol = zeros(nt,nh);
 int_S_beta = zeros(nt,nh) ; 
 
for it = 1:nt
      
 tau_mol(it,1) = 0;
 int_S_beta(it,1) = 0;
%% Calculating Pp, Pp_u, Pp_l and associated standard deviations in a loop
 for iz = 2:nh

tau_mol(it, iz) = tau_mol(it,iz-1) + 0.5 * ((alpha_m(it,iz-1) + alpha_m(it,iz)) * height_res);
int_S_beta(it,iz) =  int_S_beta(it,iz-1) + 0.5 *((beta_m(it,iz-1)*S(it,iz-1) + beta_m(it,iz)*S(it,iz)) * height_res);

Pp(it,iz) = S(it,iz) * P_sig(it,iz) * exp(-2.0*int_S_beta(it,iz)) * exp(2*tau_mol(it,iz));
sd_Pp(it,iz) = S(it,iz) * sd_P_sig(it,iz) * exp(-2.0 *int_S_beta(it,iz)) * exp(2*tau_mol(it,iz));    



Pp_single(it,iz) = S(it,iz) * P_sig_single(it,iz) * exp(-2.0*int_S_beta(it,iz)) * exp(2*tau_mol(it,iz));

sd_Pp_single(it,iz) = S(it,iz) * sd_P_sig_single(it,iz) * exp(-2.0 *int_S_beta(it,iz)) * exp(2*tau_mol(it,iz));

 end

 Pp_int{it} = Pp(it,:) ;
 Pp_int_s{it} = Pp_single(it,:);

 
% Calculate the index for io  Find the index corresponding to zo
height_o(it)= (height_norm_1(it) + height_norm_2(it))/2;
io(it)  = find(height <= height_o(it), 1, 'last' );

% Calculating slope of the signal
x = log(beta_atten(it,:));
y = height;
x_s = log(beta_atten_single(it,:));

beta_slope{it} = x;
beta_slope_s{it} = x;
height_slope{it} = y;

  
%   lm = fitlm(x,y);
alpha_o_p(it,:) = -0.5 * (diff(beta_slope{it}) ./ height_res);% +  beta_m(it,io(it)) * S(it,io(it));
% alpha_slope(it) = alpha_o_p(it,round(length(beta_slope{it})/2)); 
alpha_slope(it) = (alpha_o_p(it,io(it)))+bias;

alpha_o_p_s(it,:) = -0.5 * (diff(beta_slope_s{it}) ./ height_res);% +  beta_m(it,io(it)) * S(it,io(it));
% alpha_slope_s(it) = alpha_o_p_s(it,round(length(beta_slope_s{it})/2)); 
alpha_slope_s(it) = (alpha_o_p_s(it,io(it)))+bias;
  
% Estimating value of alpha_o 
alpha_o(it) = Ext(it,io(it)) + S(it, io(it)) .* beta_m(it, io(it));


[alphap{it}, sd_alphap{it}] = Klett_inversion(height, Pp_int{it}, sd_Pp, ...
     io(it), alpha_slope(it)) ;
[alphap_s{it}, sd_alphap_s{it}] = Klett_inversion(height, Pp_int_s{it}, sd_Pp_single, ...
     io(it), alpha_slope(it)) ; 
 

alpha(it,:) = alphap{it} - S(it,:) .* beta_m(it,:) ; % becuse there is no aerosols below the cloud in ECSIM simulations
sd_alpha(it,:) = sd_alphap{it} ;%- S(it,:) .* beta_m(it,:) ;

alpha_s(it,:) = alphap_s{it}- S(it,:) .* beta_m(it,:)   ;%-

% Calculating the corresponding beta (backscattering) values
beta(it,:) = alpha(it,:) ./ S(it,iz) ;
sd_beta(it,:) = sd_alpha(it,:) ./ S(it,iz) ;
  
clear Pp sd_Pp x y
end
       
tau_org = height_res * cumsum(Ext,1);
tau = height_res * cumsum(alpha,1);         %Optical depth
sd_tau = height_res * cumsum(sd_alpha,1) ;  %Estimated random error in optical depth
e = std(alpha);

tt=120;

figure('NumberTitle','off', ...
    'Units', 'centimeters','Position',[2 20 15 15])
p1 = plot(alpha(tt,id_cb_lidar(tt):norm_down(tt)), ...
    tau_org(tt,id_cb_lidar(tt):norm_down(tt)),'-.r',...
    Ext(tt,id_cb_lidar(tt):norm_down(tt)), ...
    tau_org(tt,id_cb_lidar(tt):norm_down(tt)),'-.b', ...
    alpha_s(tt,id_cb_lidar(tt):norm_down(tt)), ...
    tau_org(tt,id_cb_lidar(tt):norm_down(tt)), '-.k', ...
    'LineWidth',1.5,...
    'MarkerSize',10);
% ylim([0 10])
t = title('Retrieved Extinction coefficient');
set(t, 'horizontalAlignment', 'left')
set(t, 'units', 'normalized')
h1 = get(t, 'position');
set(t, 'position', [0 h1(2) h1(3)])
grid on
legend('Retrieved extinction', '"True" extinction', 'Retrieved extinction with multiple scattering correction')
xlabel('\alpha [m^{-1} sr]')
ylabel('\tau ')

%% Simple check of accuracy of alpha_o
for i = 1:nt
    accuracy(i) = alpha_slope_s(i)./Ext(i,io(i));
    bias(i) = alpha_slope_s (i) -(Ext(i,io(i)) + S(it, io(it)) .* beta_m(it, io(it)));
end
mean_bias = mean(bias)
max_accu_alpha_o = max(accuracy(~isinf(accuracy)))
min_accu_alpha_o = min(accuracy(~isinf(accuracy)))
mean_accu_alpha_o = mean(accuracy(~isinf(accuracy)))

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
 accu_alpha(it,ih) = (alpha(it,retrieved_cloud{it}(ih))./ Ext(it,retrieved_cloud{it}(ih)));
 err_alpha(it,ih) = (alpha(it,retrieved_cloud{it}(ih)) - Ext(it,retrieved_cloud{it}(ih))); 
 accu_alpha_sscat(it,ih) = (alpha_s(it,retrieved_cloud{it}(ih)) ./ Ext(it,retrieved_cloud{it}(ih)));
 err_alpha_sscat(it,ih) =(alpha_s(it,retrieved_cloud{it}(ih)) - Ext(it,retrieved_cloud{it}(ih))); 
 single_improvement(it,ih) =  (err_alpha(it,ih) - err_alpha_sscat(it,ih)) ;   
 bias_alpha(:,ih) = (sum(alpha(:,retrieved_cloud{it}(ih))) - sum(Ext(:,retrieved_cloud{it}(ih)))) ./nt ;
RMSE(:,ih) = sqrt(mean((Ext(it,retrieved_cloud{it}(ih)) - alpha(it,retrieved_cloud{it}(ih))).^2));
RMSE_ss(:,ih) = sqrt(mean((Ext(it,retrieved_cloud{it}(ih)) - alpha_s(it,retrieved_cloud{it}(ih))).^2));
end

accu_alpha_plot(it,1:retrieved_cloud{it}-1)=NaN;
accu_alpha_plot(it,max(retrieved_cloud{it})+1:end)=NaN;
   
end
accu_alpha_plot(accu_alpha_plot ==0) = NaN;
accu_alpha(accu_alpha ==0) = NaN;
err_alpha(err_alpha==0)= NaN;
accu_alpha_sscat(accu_alpha_sscat==0)= NaN;
err_alpha_sscat(err_alpha_sscat==0)= NaN;

%% Latex table
fprintf('height above cb & accu_alpha & error alpha & accu_alpha_SS & error alpha SS & improvement & bias alpha &RMSE & RMSE_ss\n')
for ii=1:6
    fprintf('%8.0f & %8.2f &  %8.6f &  %8.2f & %8.6f &  %8.6f &  %8.6f &  %8.6f&  %8.6f \\\\ \n', ...
        height(retrieved_cloud{ii}(ii)) - height(retrieved_cloud{ii}(1)) , ...
        nanmean(accu_alpha(:,ii)), ...
        nanmean(err_alpha(:,ii)), ...
        nanmean(accu_alpha_sscat(:,ii)),...
        nanmean(err_alpha_sscat(:,ii)), ...
        mean(single_improvement(:,ii)), ...
        bias_alpha(:,ii), ...
        nanmean(RMSE(ii)), ...
        nanmean(RMSE_ss(ii)))

end

%% Plot profile with horizontal errorbar based on the bias of the measurements
% figure('NumberTitle','off', ...
%     'Units', 'centimeters','Position',[2 80 15 15])
% hb1 = herrorbar(alpha_s(tt,id_cb_lidar(tt):norm_down(tt)), ...
%     tau_org(tt,id_cb_lidar(tt):norm_down(tt)), ...
%     err_alpha_sscat(tt,1:length(alpha_s(tt,id_cb_lidar(tt):norm_down(tt)))));
% set(hb1, 'LineWidth',2)
% %     bias_alpha(1:length(alpha_s(tt,id_cb_lidar(tt):norm_down(tt)))))
% hold on
% p12 = plot(Ext(tt,id_cb_lidar(tt):norm_down(tt)), ...
%     tau_org(tt,id_cb_lidar(tt):norm_down(tt)), '-k','MarkerSize',8);
% set(p12, 'LineWidth',2)
% ylim([0 80])
% xlabel('\alpha [m^{-1} sr]')
% ylabel('\tau ')
% t2 = title('Retrieved Extinction coefficient +/- bias');
% set(t2, 'horizontalAlignment', 'left')
% set(t2, 'units', 'normalized')
% h1 = get(t2, 'position');
% set(t2, 'position', [0 h1(2) h1(3)])
% grid on

% figure('NumberTitle','off', ...
%     'Units', 'centimeters','Position',[2 120 15 15])
% s = 20;
% sp1 = scatter(alpha_s(tt,id_cb_lidar(tt):norm_down(tt)),...
%     Ext(tt,id_cb_lidar(tt):norm_down(tt))',s,tau_org(tt,id_cb_lidar(tt):norm_down(tt)),'filled');
% % sp1 = plotmatrix(alpha_s(:,id_cb_lidar:norm_down),...
% %     Ext(:,id_cb_lidar:norm_down));
% ylim([0 0.05])
% xlim([0 0.05])
% xlabel('\alpha retrieved [m^{-1} sr]')
% ylabel('\alpha modelled [m^{-1} sr] ')
% t3 = title('Retrieved extinction vs Modelled extinction');
% refline(1,0) 
% lsline
% colorbar
% set(t3, 'units', 'normalized')

% %% Create scatter plots
% for it = 1:nt
%     for ir = 1:length(retrieved_cloud{it})
%     data_alpha{ir} = alpha_s( :, retrieved_cloud{it}(ir)) ; 
%     data_ext{ir} = Ext(:,retrieved_cloud{it}(ir));
%     data_tau{ir} = tau_org(:,retrieved_cloud{it}(ir));
%     data_height{ir} = height(retrieved_cloud{it}(ir)) - height(retrieved_cloud{it}(1));
%     end
% end
% 
% %% Plot scatter plot of retrieved and modelled extionction coefficient
% TitleFigure=['scatter plot of retrieved and modelled extionction coefficient'];
% figure('name', TitleFigure, 'NumberTitle','off', ...
%     'Units','centimeters','Position',[10 30 15 15]);
% s = 25;
% for is = 1:6
%  cx(is) = subplot(2,3,is) ; 
%  scat_fig(is) = scatter(data_alpha{is},data_ext{is}, s,  data_tau{is} , 'fill');
%  hold on
%  refline(1,0) 
% %  set(gca, 'YScale', 'log')
% %  set(gca, 'XScale', 'log')
% %  set(gca, 'XTickLabel', [0.5 1 1.5 2 2.5 3]); 
% %  set(gca, 'FontSize',8)
%  caxis([min(min([data_tau{:}])) max(max([data_tau{:}]))])
%  xlim([0 0.05])
%  ylim([0 0.05])
% %  a=get(gca,'XLim');
% %  x=max(a)-(max(a)-1.1*min(a));
% %  b=get(gca,'YLim');
% %  y=max(b)-(max(b)-min(b))/10;
% %  text(x,y,[...%'\itr = ',num2str(correlations_reff_beta{i}(2),'%.2f'),...
% %     ...%'\newline','\itr^2 = ',num2str(mdl_lidar_reff{i}.Rsquared.Ordinary,'%.2f'),...
% %     ...%'\newline','\itm = ',num2str(poly_reff_beta{i}(1),'%.2f') ...
% %     '\newline','\itACI_{r} = ',num2str((poly_reff_aci{i}(1))*-1,'%.2f')],'FontSize',10);
%  title(['Height above cloud base ' num2str(data_height{is})],'FontSize',10,'FontWeight','normal')
% end
% H=labelEdgeSubPlots('\alpha modelled [m^{-1} sr]','\alpha retrieved [m^{-1} sr]');
% h2=colorbar;
% set(h2, 'Position', [.925 .11 .0181 .81])
% h2.Label.String = '\tau';
% h2.Label.FontSize = 8;
% h2.Label.Rotation = 0; 
% h2.Label.Position = [2 15 0];
figure('NumberTitle','off', ...
    'Units', 'centimeters','Position',[2 20 15 15])
p2 = plot((P_sig(tt,id_cb_lidar(tt):norm_down(tt))./P_sig_single(tt,id_cb_lidar(tt):norm_down(tt))), ...
    tau_org(tt,id_cb_lidar(tt):norm_down(tt)),'-.r','LineWidth',1.5, ...
    'MarkerSize',10);
% ylim([0 10])
t = title('P_{sig}/P_{sig_{single}}');
set(t, 'horizontalAlignment', 'left')
set(t, 'units', 'normalized')
h1 = get(t, 'position');
set(t, 'position', [0 h1(2) h1(3)])
grid on
% legend('Retrieved extinction', '"True" extinction', 'Retrieved extinction with multiple scattering correction')
xlabel(' ')
ylabel('\tau ')

figure('NumberTitle','off', ...
    'Units', 'centimeters','Position',[2 20 15 15])
p3 = plot(depol_ratio(tt,id_cb_lidar(tt):norm_down(tt)), ...
    tau_org(tt,id_cb_lidar(tt):norm_down(tt)),'-.r','LineWidth',1.5, ...
    'MarkerSize',10);
% ylim([0 10])
t = title('Depolarisation ration');
set(t, 'horizontalAlignment', 'left')
set(t, 'units', 'normalized')
h1 = get(t, 'position');
set(t, 'position', [0 h1(2) h1(3)])
grid on
% legend('Retrieved extinction', '"True" extinction', 'Retrieved extinction with multiple scattering correction')
xlabel('depolarisation ratio')
ylabel('\tau ')