% Lidar inversion - module for simulated data
clear all 
format long

% % 
% load('les_1-merged_with_driz_5m.mat')
% load('true_les_1-merged_with_driz_5m.mat')
load('fire_27_merged_with_mask_v6.mat')
load('fire_27_merged_with_mask_v6_ss.mat')
load('true_fire_27_merged_with_mask_v6.mat')


beta_atten = lid_total ; % Original unit is 1/m/sr, change to 1/km/sr
% s = RandStream('mt19937ar','seed',5489);
% beta2 = awgn(beta_atten, 5, 0, s);
% beta_atten = beta2;
sd_beta_atten = 0.2 .* beta_atten ;

 nh = nz;   
 
 nt = nx;
 time = linspace(1,nt,nt);
 height = height(1,:) .* 1.0e3;
 height_res = height(2) - height(1);  %Height resolution in km


% Conver Range Corrected Signal to Power
P_sig = zeros(nt,nh);

 for it = 1:nt
 P_sig(it,:) = beta_atten(it,:) ./ height.^2;
 end
 sd_P_sig = 0.2 .* P_sig;
% s = RandStream('mt19937ar','seed',5489);
% P_noise = awgn(P_sig, 3, 0, s);
% 
% for it = 1:nt
%     
%  beta_atten(it,:) = P_noise(it,:) .* height.^2;
% end
 
 
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

    height_norm_1(it) = height(norm_down(it));
    height_norm_2(it) = height(norm_up(it)) ;

    height_o(it) = (height_norm_1(it) + height_norm_2(it))/2;
    io(it)  = find(height <= height_o(it), 1, 'last' );

 end

%% Corecting signal for the multiple scattering
perp_int = zeros(nt,nh);
para_int = zeros(nt,nh);
perp_sum = zeros(nt,nh);
para_sum = zeros(nt,nh);

int_MS = beta_atten;
beta_atten_single = beta_atten;
P_sig_single = P_sig ;

corrected_signal = beta_atten;
for it= 1:nt

    perp_sum (it,id_cb_lidar(it)-1:id_ct_radar(it)) = cumsum( perp(it, ...
        id_cb_lidar(it)-1:id_ct_radar(it)),2);
    para_sum (it,id_cb_lidar(it)-1:id_ct_radar(it)) = cumsum( para(it, ...
        id_cb_lidar(it)-1:id_ct_radar(it)),2);
    
    for jz = 1:length(id_cb_lidar(it)-1:id_ct_radar(it))
    perp_int(it,id_cb_lidar(it)-1+jz) = perp_sum (it,id_cb_lidar(it)-1+jz)* height_res;
   
    para_int(it,id_cb_lidar(it)-1+jz) = para_sum (it,id_cb_lidar(it)-1+jz)* height_res;
    
    int_MS(it,id_cb_lidar(it)-1+jz) = (perp_sum(it,id_cb_lidar(it)-1+jz) + ...
        para_sum(it,id_cb_lidar(it)-1+jz)) * height_res ;
    end 
        
end
depol_ratio = perp_int ./ para_int ;
correction_factor = (1 - depol_ratio) ./ (1 + depol_ratio);
correction_factor(isnan(correction_factor)) = 1 ;
int_SS = int_MS .* correction_factor ;

for it= 1:nt
   for jjz = 1:length(id_cb_lidar(it)-1:id_ct_radar(it))
    beta_atten_single(it,id_cb_lidar(it)-1+jjz) = (int_SS(it,id_cb_lidar(it)-1+jjz) - ...
        int_SS(it,id_cb_lidar(it)-1+jjz-1)) ./ height_res ; 
    P_sig_single(it,id_cb_lidar(it)-1+jjz) = beta_atten_single(it,id_cb_lidar(it)-1+jjz) ...
        ./ (height(id_cb_lidar(it)-1+jjz).^2) ;
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


 
% Calculate the index for io  Find the index corresponding to zo
height_o(it)= (height_norm_1(it) + height_norm_2(it))/2;
io(it)  = find(height <= height_o(it), 1, 'last' );

% Calculating slope of the signal
x = log(beta_atten(it,(norm_down:norm_up)));
y = height(norm_down:norm_up);
x_s = log(beta_atten_single(it,(norm_down:norm_up)));

beta_slope{it} = x;
beta_slope_s{it} = x;
height_slope{it} = y;

  
%   lm = fitlm(x,y);
alpha_o_p(it,:) = -0.5 * (diff(beta_slope{it}) ./ height_res);% +  beta_m(it,io(it)) * S(it,io(it));
% alpha_slope(it) = alpha_o_p(it,round(length(beta_slope{it})/2)); 
alpha_slope(it) = mean(alpha_o_p(it,:));

alpha_o_p_s(it,:) = -0.5 * (diff(beta_slope_s{it}) ./ height_res);% +  beta_m(it,io(it)) * S(it,io(it));
% alpha_slope_s(it) = alpha_o_p_s(it,round(length(beta_slope_s{it})/2)); 
alpha_slope_s(it) = mean(alpha_o_p_s(it,:));
  
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

tt=35;

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
% hold on
% p2 = plot(Ext(tt,id_cb_lidar(tt):norm_down(tt)), ...
%     tau_org(tt,id_cb_lidar(tt):norm_down(tt)), ...
%     '-.b','MarkerSize',8);
% hold on 
% p3 = plot(alpha_s(tt,id_cb_lidar(tt):norm_down(tt)), ...
%     tau_org(tt,id_cb_lidar(tt):norm_down(tt)), ...
%     '-.k','MarkerSize',8);
% set(p1, 'LineWitdth' ,2 )
% hold on 
ylim([0 80])
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
    accuracy(i) = alpha_slope(i)./Ext(i,io(i));
end
max(accuracy(~isinf(accuracy)))
min(accuracy(~isinf(accuracy)))
mean(accuracy(~isinf(accuracy)))

%% Calculate statistical data and make tables

for it = 1:nt
retrieved_cloud{it} = id_cb_lidar(it):norm_down(it);
height_cloud{it} = height(retrieved_cloud{it});
end


for it = 1:nt
for ih = 1:length(retrieved_cloud{it})
 accu_alpha(it,ih) = (alpha(it,retrieved_cloud{it}(ih))./ Ext(it,retrieved_cloud{it}(ih)));
 err_alpha(it,ih) = (alpha(it,retrieved_cloud{it}(ih)) - Ext(it,retrieved_cloud{it}(ih))); 
 accu_alpha_sscat(it,ih) = (alpha_s(it,retrieved_cloud{it}(ih)) ./ Ext(it,retrieved_cloud{it}(ih)));
 err_alpha_sscat(it,ih) =(alpha_s(it,retrieved_cloud{it}(ih)) - Ext(it,retrieved_cloud{it}(ih))); 
 single_improvement(it,ih) =  (err_alpha(it,ih) - err_alpha_sscat(it,ih)) ;   
 bias_alpha(:,ih) = (sum(alpha(:,retrieved_cloud{it}(ih))) - sum(Ext(:,retrieved_cloud{it}(ih)))) ./nt ;
end
   
end

%% Latex table
fprintf('height above cb & accu_alpha & error alpha & accu_alpha_SS & error alpha SS & improvement & bias alpha \n')
for ii=1:6
    fprintf('%8.0f & %8.6f &  %8.6f &  %8.6f & %8.6f &  %8.6f &  %8.6f  \\\\ \n', ...
        height(retrieved_cloud{ii}(ii)) - height(retrieved_cloud{ii}(1)) , ...
        mean(accu_alpha(:,ii)), ...
        mean(err_alpha(:,ii)), ...
        mean(accu_alpha_sscat(:,ii)),...
        mean(err_alpha_sscat(:,ii)), ...
        mean(single_improvement(:,ii)), ...
        bias_alpha(:,ii))

end

%% Plot profile with horizontal errorbar based on the bias of the measurements
figure('NumberTitle','off', ...
    'Units', 'centimeters','Position',[2 80 15 15])
hb1 = herrorbar(alpha_s(tt,id_cb_lidar(tt):norm_down(tt)), ...
    tau_org(tt,id_cb_lidar(tt):norm_down(tt)), ...
    err_alpha_sscat(tt,1:length(alpha_s(tt,id_cb_lidar(tt):norm_down(tt)))));
set(hb1, 'LineWidth',2)
%     bias_alpha(1:length(alpha_s(tt,id_cb_lidar(tt):norm_down(tt)))))
hold on
p12 = plot(Ext(tt,id_cb_lidar(tt):norm_down(tt)), ...
    tau_org(tt,id_cb_lidar(tt):norm_down(tt)), '-k','MarkerSize',8);
set(p12, 'LineWidth',2)
ylim([0 80])
xlabel('\alpha [m^{-1} sr]')
ylabel('\tau ')
t2 = title('Retrieved Extinction coefficient +/- bias');
set(t2, 'horizontalAlignment', 'left')
set(t2, 'units', 'normalized')
h1 = get(t2, 'position');
set(t2, 'position', [0 h1(2) h1(3)])
grid on

figure('NumberTitle','off', ...
    'Units', 'centimeters','Position',[2 120 15 15])
sp1 = scatter(alpha_s(tt,id_cb_lidar(tt):norm_down(tt)),...
    Ext(tt,id_cb_lidar(tt):norm_down(tt))');
