format long
close all
% Load data from ECSIM simulation
load('fire_27_merged_with_mask_v7.mat')
load('fire_27_merged_with_mask_v7_ss.mat')
load('true_fire_27_merged_with_mask_v7.mat')
load('correction_factor.mat')

%% Select a specific part of the simulation
% Simulation runs from 0 to 250
start = 1;
stop = 250;
% Select needed varbles in the selected period
Ba_MS = lid_total(start:stop,1:end-1); % Original unit is 1/m/sr
depol_MS = lin_depol(start:stop,1:end-1);

Ba_SS = lid_total_ss(start:stop,2:end);
depol_SS = lin_depol_ss(start:stop,2:end);

% Divide multiple signal into parallel and perpendicular
para_ms = Ba_MS ./ (depol_MS + 1) ;
perp_ms = Ba_MS - para_ms;

Ze = Ze(start:stop,1:end-1);
Ext = Ext(start:stop,2:end);
P = P(start:stop,1:end-1);
T = T(start:stop,1:end-1);
LWP = LWP(start:stop,1:end-1);

sd_beta_atten = 0.2 .* Ba_MS ;

nz = size(lid_total(start:stop,1:end-1),2);   
nt = size(lid_total(start:stop,1:end-1),1);
along_track = linspace(1,nt,nt);
z = height(1,1:end-1) .* 1.0e3;
res = z(2) - z(1);  %Height resolution in meters

%% Estimating cloud base height and defining normalisation interval
 [Peak_P, I_max] = max(Ba_MS, [], 2) ;
 [Peak_perp_ms, perp_max] = max(perp_ms, [], 2) ;
 
 id_cb_lidar = zeros(1,nt); 
 id_ct_radar = zeros(1,nt); 
 norm_down = zeros(1,nt);
 norm_up = zeros(1,nt);
 z_norm_1 = zeros(1,nt);
 z_norm_2 = zeros(1,nt);
 z_o = zeros(1,nt);
 io = zeros(1,nt);

 for it = 1:nt
     
    id_cb_lidar(it) = find (perp_ms(it,:) >= (Peak_perp_ms(it)/10), 1,'first'); 
    id_ct_radar(it) = find((Ze(it,:)) >= 1.e-3, 1, 'last');

    norm_down(it) = I_max(it)+4;
    norm_up(it) = norm_down(it) +2;    

    z_norm_1(it) = z(norm_down(it));
    z_norm_2(it) = z(norm_up(it)) ;

    z_o(it) = (z_norm_1(it) + z_norm_2(it))/2;
    io(it)  = find(z <= z_o(it), 1, 'last' );

 end

%% Correct signal from the multiple scattering simulation
perp_int = zeros(nt,nz);
para_int = zeros(nt,nz);
int_MS = zeros(nt,nz);

Ba_SS_corr = Ba_MS .* res;
ms_down = id_cb_lidar;
ms_up = id_ct_radar ;

for it= 1:nt
       
    perp_int(it,ms_down(it):ms_up(it))= cumsum(perp_ms(it,...
        ms_down(it):ms_up(it)),2).* res;
    para_int(it,ms_down(it):ms_up(it))= cumsum(para_ms(it,...
        ms_down(it):ms_up(it)),2).* res;
%     perp_int(it,ms_down(it)) = perp_ms(it,ms_down(it)) .* res;
%     para_int(it,ms_down(it)) = para_ms(it,ms_down(it)) .* res; 
%     int_MS(it, ms_down(it)) = perp_int(it,ms_down(it)) + ...
%                               para_int(it,ms_down(it));        
%     
%     for jz = ms_down(it)+1:ms_up(it)
%         perp_int(it,jz) = perp_int(it,jz-1) + 0.5.*(perp_ms(it,jz-1).*res ...
%             + perp_ms(it,jz).* res);
%         para_int(it,jz) = para_int(it,jz-1) + 0.5.*(para_ms(it,jz-1).*res ...
%             + para_ms(it,jz).* res); 
%         int_MS(it,jz) = perp_int(it,jz) + para_int(it,jz);
%     end    
        
    int_MS(it,ms_down(it):ms_up(it)) = ...
        (perp_int(it,ms_down(it):ms_up(it)) + ...
        para_int(it,ms_down(it):ms_up(it)));  
end
depol_ratio = perp_int ./ para_int ;
correction_factor = ((1 - depol_ratio) ./ (1 + depol_ratio)).^2;
correction_factor(isnan(correction_factor)) = 1 ;
f = 0;
A_p = (1-correction_factor).*f + correction_factor;
% int_SS = int_MS.*AS(:,1:end-1);
int_SS = int_MS.*A_p;

dM = correction_factor;
dAS = AS(:,1:end-1);


for it= 1:nt
   
% Ba_SS_corr(it,ms_down(it):ms_up(it)) = diff(int_SS(it,ms_down(it)-1:ms_up(it))) ./ res;
    
  Ba_SS_corr(it,ms_down(it))= int_SS(it,ms_down(it));
 for jjz = ms_down(it)+1:ms_up(it)
  Ba_SS_corr(it,jjz) =  (int_SS(it,jjz+1)) - (int_SS(it,jjz-1));
 end
%  dM(it,1) = correction_factor(it,1)
 for j = 2:length(z)-1
     dM(it,j) = dM(it,j+1) - dM(it,j);
     dAS(it,j) = dAS(it,j+1) - dAS(it,j);
 end
end
dMdz = dM ./ res;
dASdz = dAS ./ res;

figure
plot(dMdz(160,1:40),z(1:40), 'b', ...
    dASdz(160,1:40),z(1:40), 'k')

Ba_SS_corr = Ba_SS_corr./res;

Ba_SS_corr2 = correction_factor .* Ba_MS + int_MS .* dMdz;
figure
plot(Ba_SS_corr(120,:),z, '-k', ...
    Ba_SS(120,:),z,'-.r', ...
    Ba_SS_corr2(120,:),z, '-.b', ...
    Ba_MS(120,:),z, '.-m')
ylim([0 600])


figure 
plot(correction_factor(100,:)./AS(100,1:end-1),z,'.-b')
ylim([0 600])
hold on
 line('XData', [0 1], 'YData', [z(ms_down(100)) z(ms_down(100))], ...
    'LineWidth', 1.0, 'LineStyle', '-.', 'Color', 'm')

%% Find Rayleigh extinction and backscattering cross sections 
%  This method is based on D.R. Bates, 1984 and is valid for wavelengths
%  from about 300 to 1000 nm.

% Calculate the rayleigh extinction coefficient in cm^2/sr-1. 
% Wavelength is in nm. Must multiply by 1.e+5 if dealing with km.
 wavelen = 350;

    work = wavelen*1.0e-3 ;
     if wavelen < 550
            x = 0.389 * work + 0.09426/work - 0.3228 ;
      else  x = 0.04 ;
     end
    sigma = (4.02e-28/ (work ^(4 + x))); 

    % Calculate King's Factor
    FK = 1.225e-13*wavelen^4 - 3.911e-10*wavelen^3 + 4.6100e-7*wavelen^2 ...
        -2.410e-4*wavelen + 1.095 ;  
    e=(FK-1.0)*9.0/2.0;
    beta = 3.0 / (4.0*pi) * (180.0+28.0*e) / (360.0+82*e) * sigma;
    % Rayleigh extinction cross section
    sigma = sigma*1.e+5;
    % Rayleigh backscattering cross section
    beta = beta*1.e+5;    

%% Calculate Rayleigh extinction and backscatter
% Calculate atmospheric density 
k_boltz = 1.380658e-23; % Nm/K
rho_atm = P./T./(k_boltz.*1.e+4); 
% 1.e+4 is to compensate for m --> cm and mb --> pa
beta_ray = rho_atm.* beta .*1.e-3; %1/km
ext_ray = rho_atm.*sigma .*1.e-3; %1/km

%% Calculate optical thickness and Rayleigh optical thickness
extp=Ext;
% 
ext_rayp=ext_ray;
beta_rayp=beta_ray;
 
% Calculate optical thickness
  tau=zeros(nt,nz);
  tau_ray=zeros(nt,nz);
  int_beta_ray=zeros(nt,nz);
    
  for it = 1:nt
    tau(it,1) =  res.*(extp(it,1)+ext_rayp(it,1));
    int_beta_ray(it,1)= res.*beta_rayp(it,1);
    tau_ray(it,1) = res.*ext_rayp(it,1);
      for j = 2:nz
        tau(it,j) = tau(it,j-1)+ 0.5.*res.* ...
         (extp(it,j-1)+ext_rayp(it,j-1)+extp(it,j)+ext_rayp(it,j));
        int_beta_ray(it,j) = int_beta_ray(it,j-1)+...
         res*(beta_rayp(it,j-1)+beta_rayp(it,j))*0.5;
        tau_ray(it,j)=tau_ray(it,j-1)+res*(ext_rayp(it,j-1)+...
            ext_rayp(it,j))*0.5;
      end     
  end
  
%% Add Gaussian error to th signal


  
%% Calucalte single and multiple attenuated backscatter 
%  corrected for the Reayleigh   
  S=16.0;
  Ba_msp=S.*Ba_MS.*exp(2.*tau_ray).*exp(-2.*S.*int_beta_ray);
  Ba_ssp=S.*Ba_SS.*exp(2.*tau_ray).*exp(-2.*S.*int_beta_ray);
  Ba_ssp_corr = S.*Ba_SS_corr2.*exp(2.*tau_ray).*exp(-2.*S.*int_beta_ray);
  
%   figure
% plot(Ba_ssp_corr(45,:),z, 'k', ...
%     Ba_ssp(45,:),z,'r')
% ylim([0 600])
    
%% Caluclate resolution correction factor
%  Here calculated from 'real' extinction from ECSIM

rfac=1.0./(2.0.*res.*(extp+ext_rayp)).*(1.0-exp(-2.0.*(extp+ext_rayp).*res));


%% Perform signal inversion
alpha_ret_ss = zeros(nt,nz);
alpha_ret_ss_nc = zeros(nt,nz);
alpha_ret_ms = zeros(nt,nz);
alpha_ret_ms_nc = zeros(nt,nz);
alpha_ret_ms2ss = zeros(nt,nz);

 for it = 1:nt
[alpha_ss_slope(it,:)] = slope(z,Ba_ssp(it,:));
[alpha_ms_slope(it,:)] = slope(z,Ba_msp(it,:));
[alpha_ss_slope_corr(it,:)] = slope(z,Ba_ssp_corr(it,:));


[alpha_ret_ss_nc(it,:)] = simple_inversion(io(it),...
    (extp(it,:)+S.*beta_ray(it,:)),z, Ba_ssp(it,:),0,0);

[alpha_ret_ss(it,:)] = simple_inversion(io(it),...
    alpha_ret_ss_nc(it,:),z, Ba_ssp(it,:),1,1);

[alpha_ret_ms_nc(it,:)] = simple_inversion(io(it),...
     (alpha_ms_slope(it,:)),z, Ba_msp(it,:),0,0);

[alpha_ret_ms(it,:)] = simple_inversion(io(it),...
    alpha_ret_ms_nc(it,:),z, Ba_msp(it,:),1,1);

[alpha_ret_ms2ss_nc(it,:)] = simple_inversion(io(it),...
    alpha_ss_slope_corr(it,:),z, Ba_ssp_corr(it,:),0,0);

[alpha_ret_ms2ss(it,:)] = simple_inversion(io(it),...
    alpha_ret_ms2ss_nc(it,:),z, Ba_ssp_corr(it,:),1,1);

alpha_ret_ss(it,:) = alpha_ret_ss(it,:)-S.*beta_ray(it,:);
alpha_ret_ss_nc(it,:) = alpha_ret_ss_nc(it,:)-S.*beta_ray(it,:);
alpha_ret_ms(it,:) = alpha_ret_ms(it,:) -S.*beta_ray(it,:);
alpha_ret_ms_nc(it,:) = alpha_ret_ms_nc(it,:) -S.*beta_ray(it,:);
alpha_ret_ms2ss(it,:) = alpha_ret_ms2ss(it,:) -S.*beta_ray(it,:);
 end
 
figure
tt=205;
 plot(extp(tt,:),z,'.-k', ...
     alpha_ret_ss(tt,:),z, '.-g', ...
     alpha_ret_ss_nc(tt,:),z,'.-r', ...
     alpha_ret_ms(tt,:),z,'.-b', ...
     alpha_ret_ms_nc(tt,:),z,'.-y', ...
     alpha_ret_ms2ss(tt,:),z,'.-c')
 line('XData', [0 0.05], 'YData', [z(io(tt)) z(io(tt))], ...
    'LineWidth', 1.0, 'LineStyle', '-.', 'Color', 'm')
 xlim([0 0.04])
 ylim([0 500])
 l = legend('True extinction', ...
    'Retrieved extinction SS with RES correction', ...
    'Retrieved extinction SS without RES correction', ...
    'Retrieved extinction MS with RES correction', ...
    'Retrieved extinction MS without RES correction', ...
    'Retrieved extinction with SS & RES correction ', ...
    'Normalisation height');

figure
plot(correction_factor(tt,1:20),z(1:20),'.-b', ...
    AS(tt,1:20),z(1:20),'.-r',...
    A_p(tt,1:20),z(1:20),'.-g')
 line('XData', [0.6 1], 'YData', [z(io(tt)) z(io(tt))], ...
    'LineWidth', 1.0, 'LineStyle', '-.', 'Color', 'm')
 line('XData', [0.6 1], 'YData', [z(I_max(tt)) z(I_max(tt))], ...
    'LineWidth', 1.0, 'LineStyle', '-.', 'Color', 'y')

% figure
%  plot(extp(tt,:),tau(tt,:),'.-k', ...
%      alpha_ret_ss(tt,:),tau(tt,:), '.-g', ...
%      alpha_ret_ss_nc(tt,:),tau(tt,:),'.-r')
%  line('XData', [0 0.05], 'YData', [z(norm_down(tt)) z(norm_down(tt))], ...
%     'LineWidth', 1.0, 'LineStyle', '-.', 'Color', 'b')
%  xlim([0 0.04])
%  ylim([0 1])
%  l = legend('True extinction', ...
%     'Retrieved extinction with RES correction', ...
%     'Retrieved extinction without RES correction', ...
%     'Begining of the normalisation invterval');

% figure
% p1 = pcolor(along_track,z,Ba_MS');
% set(p1, 'EdgeColor', 'none');
% set(gca,'Ydir','normal');
% set(gca,'ylim',[0 500])
% ylabel('Height [m]')
% hold on
% l2 = plot(along_track,z(norm_down),'k','LineWidth',1.5);
% % hold on
% set(gca, 'FontSize',10)
% c=colorbar('southoutside');
% c.Label.String = 'ATB [1/m/sr]';
% c.Label.FontSize = 10;
% t = title('Attenuated backscatter', 'FontSize',12,'FontWeight','normal');
% set(t, 'horizontalAlignment', 'left')
% set(t, 'units', 'normalized')
% h1 = get(t, 'position');
% set(t, 'position', [0 h1(2) h1(3)])
% grid on
% h_legend=legend([l2,l3], 'Start of the normalisation interval', 'Cloud base estimate');
% set(h_legend,'FontSize',10);
%% Error calculation
%% Simple check of accuracy of alpha_o
for i = 1:nt
    accuracy_ss(i) = alpha_ss_slope(i,io(i))./Ext(i,io(i));
    accuracy_ms(i) = alpha_ms_slope(i,io(i))./Ext(i,io(i));
    bias(i) = alpha_ss_slope (i) - (Ext(i,io(i)));
end
mean_bias = mean(bias)
max_accu_alpha_o = max(accuracy_ss(~isinf(accuracy_ss)))
min_accu_alpha_o = min(accuracy_ss(~isinf(accuracy_ss)))
mean_accu_alpha_o = mean(accuracy_ss(~isinf(accuracy_ss)))

max_accu_alpha_o_ms = max(accuracy_ms(~isinf(accuracy_ms)))
min_accu_alpha_o_ms = min(accuracy_ms(~isinf(accuracy_ms)))
mean_accu_alpha_o_ms = mean(accuracy_ms(~isinf(accuracy_ms)))
% 
%% Calculate statistical data and make tables

for it = 1:nt
retrieved_cloud{it} = id_cb_lidar(it):norm_down(it);
height_cloud{it} = z(retrieved_cloud{it});
end


for it = 1:nt
    for iih = 1:nz
      accu_alpha_plot(it,iih) = (alpha_ret_ms2ss(it,iih) ./ Ext(it,iih));
      err_alpha_MS_corr_plot(it,iih) = (abs(alpha_ret_ms2ss(it,iih) ...
     - Ext(it,iih))./ Ext(it,iih))*100; 
      err_alpha_MS_plot(it,iih) = (abs(alpha_ret_ms2ss(it,iih) ...
     - Ext(it,iih))./ Ext(it,iih))*100; 
     end  
    
for ih = 1:length(retrieved_cloud{it})
 accu_alpha_ss(it,ih) = (alpha_ret_ss(it,retrieved_cloud{it}(ih))...
     ./ Ext(it,retrieved_cloud{it}(ih))) *100;
 err_alpha_ss(it,ih) = (abs(alpha_ret_ss(it,retrieved_cloud{it}(ih)) ...
     - Ext(it,retrieved_cloud{it}(ih)))./ Ext(it,retrieved_cloud{it}(ih)))*100;
 
  accu_alpha_ms(it,ih) = (alpha_ret_ms(it,retrieved_cloud{it}(ih))...
     ./ Ext(it,retrieved_cloud{it}(ih))) *100;
  err_alpha_ms(it,ih) = (abs(alpha_ret_ms(it,retrieved_cloud{it}(ih)) ...
     - Ext(it,retrieved_cloud{it}(ih)))./ Ext(it,retrieved_cloud{it}(ih)))*100;
 
 accu_alpha_ms_ss(it,ih) = (alpha_ret_ms2ss(it,retrieved_cloud{it}(ih))...
     ./ Ext(it,retrieved_cloud{it}(ih))) *100;
 err_alpha_ms_ss(it,ih) = (abs(alpha_ret_ms2ss(it,retrieved_cloud{it}(ih)) ...
     - Ext(it,retrieved_cloud{it}(ih)))./ Ext(it,retrieved_cloud{it}(ih)))*100;
 
 
%  accu_alpha_sscat(it,ih) = (alpha_(it,retrieved_cloud{it}(ih)) ...
%      ./ Ext(it,retrieved_cloud{it}(ih)))*100;
%  err_alpha_sscat(it,ih) =(abs(alpha_SS(it,retrieved_cloud{it}(ih)) ...
%      - Ext(it,retrieved_cloud{it}(ih)))./ Ext(it,retrieved_cloud{it}(ih)))*100; 
%  
% accu_alpha_c(it,ih) = (alpha_corr(it,retrieved_cloud{it}(ih)) ...
%      ./ Ext(it,retrieved_cloud{it}(ih)))*100;
%  err_alpha_c(it,ih) = (abs(alpha_corr(it,retrieved_cloud{it}(ih)) ...
%      - Ext(it,retrieved_cloud{it}(ih)))./ Ext(it,retrieved_cloud{it}(ih)))*100; 
%  
% 
%  accu_alpha_sscat_c(it,ih) = (alpha_SS_corr(it,retrieved_cloud{it}(ih)) ...
%      ./ Ext(it,retrieved_cloud{it}(ih)))*100;
%  err_alpha_sscat_c(it,ih) = (abs(alpha_SS_corr(it,retrieved_cloud{it}(ih)) ...
%      - Ext(it,retrieved_cloud{it}(ih)))./ Ext(it,retrieved_cloud{it}(ih)))*100; 
 

end

accu_alpha_plot(it,1:retrieved_cloud{it}-1)=NaN;
accu_alpha_plot(it,max(retrieved_cloud{it}):end)=NaN;
err_alpha_MS_corr_plot(it,1:retrieved_cloud{it}-1)=NaN;
err_alpha_MS_corr_plot(it,max(retrieved_cloud{it}):end)=NaN;

err_alpha_MS_plot(it,1:retrieved_cloud{it}-1)=NaN;
err_alpha_MS_plot(it,max(retrieved_cloud{it}):end)=NaN;
   
end


accu_alpha_ss(accu_alpha_ss ==0) = NaN;
err_alpha_ss(err_alpha_ss==0)= NaN;
accu_alpha_ms( accu_alpha_ms ==0) = NaN;
err_alpha_ms(err_alpha_ms==0)= NaN;
accu_alpha_ms_ss(accu_alpha_ms_ss==0) = NaN;
 err_alpha_ms_ss( err_alpha_ms_ss==0)= NaN;


%% Latex table
fprintf('m from cb & ACU   & ERR  -- ECSIM SS retrieved \\ \n')
for ii=1:length(retrieved_cloud{1})
    fprintf('%8.0f & %3.2f%% &  %3.2f%%  \\\\ \n', ...
        z(retrieved_cloud{1}(ii)) - z(retrieved_cloud{1}(1)) , ...
        nanmean(accu_alpha_ss(:,ii)), ...
        nanmean(err_alpha_ss(:,ii)))
end

fprintf('m from cb & ACU   & ERR  -- ECSIM MS with RES correction \\ \n')
for ii=1:length(retrieved_cloud{1})
    fprintf('%8.0f & %3.2f%% &  %3.2f%%  \\\\ \n', ...
        z(retrieved_cloud{1}(ii)) - z(retrieved_cloud{1}(1)) , ...
        nanmean(accu_alpha_ms(:,ii)), ...
        nanmean(err_alpha_ms(:,ii)))
end

fprintf('m from cb & ACU   & ERR  -- ECSIM MS with SS & RES correction \\ \n')
for ii=1:length(retrieved_cloud{1})
    fprintf('%8.0f & %3.2f%% &  %3.2f%%  \\\\ \n', ...
        z(retrieved_cloud{1}(ii)) - z(retrieved_cloud{1}(1)) , ...
        nanmean(accu_alpha_ms_ss(:,ii)), ...
        nanmean(err_alpha_ms_ss(:,ii)))
end
 