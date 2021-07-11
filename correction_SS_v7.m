%% Single scattering correction

load('fire_27_merged_with_mask_v7.mat')
load('fire_27_merged_with_mask_v7_ss.mat')

% Select needed varbles in the selected period
% beta_atten = lid_total(:,1:end-1);% Original unit is 1/m/sr, change to 1/km/sr
% lin_depol = lin_depol(:, 1:end-1);
% beta_atten_ss = lid_total_ss(:,2:end);
% lin_depol_ss = lin_depol_ss(:,2:end);

beta_atten = lid_total;% Original unit is 1/m/sr, change to 1/km/sr
% lin_depol = lin_depol;
beta_atten_ss = lid_total_ss;
% lin_depol_ss = lin_depol_ss;

sd_beta_atten = 0.2 .* beta_atten ;
sd_beta_atten_ss = 0.2 .* beta_atten_ss; 

height = height(1,1:end-1) .* 1.0e3;
height_ss = height_ss(1,2:end) .* 1.0e3;

nh = size(lid_total,2);   
nt = size(lid_total,1);
time = linspace(1,nt,nt);
 
height_res = height(2) - height(1);  %Height resolution in meters

%% Divide signal into parallel and perpendicular
para = beta_atten ./ (lin_depol + 1) ;
perp = beta_atten - para;

para_ss = beta_atten_ss ./ (lin_depol_ss +1);
perp_ss = beta_atten_ss - para_ss;

[Peak_P, I_max] = max(beta_atten, [], 2) ;
[Peak_perp, perp_max] = max(perp, [], 2) ;
%% Estimating cloud base height and defining normalisation interval
 for it = 1:nt
    id_cb_lidar(it) = find (perp(it,:) >= (Peak_perp(it)/10), 1,'first'); 
    norm_down(it) = I_max(it)+16;
    norm_up(it) = norm_down(it) +2 ;    
 end   

%% Corecting signal for the multiple scattering

perp_sum = zeros(size(perp));
para_sum = zeros(size(para));
int_MS = zeros(size(beta_atten));

perp_sum_ss = zeros(size(perp_ss));
para_sum_ss = zeros(size(para_ss));
int_SS = zeros(size(beta_atten_ss));

ms_down = id_cb_lidar;
ms_up = norm_up;

for it= 1:nt
    
    perp_sum (it,ms_down(it):ms_up(it)) = cumsum( perp(it, ...
        ms_down(it):ms_up(it)),2);
    para_sum (it,ms_down(it):ms_up(it)) = cumsum( para(it, ...
        ms_down(it):ms_up(it)),2);
    
    perp_sum_ss (it,ms_down(it):ms_up(it)) = cumsum( perp_ss(it, ...
        ms_down(it):ms_up(it)),2);
    para_sum_ss (it,ms_down(it):ms_up(it)) = cumsum( para_ss(it, ...
        ms_down(it):ms_up(it)),2);
    
    perp_int(it,ms_down(it))= perp_sum(it,ms_down(it))* height_res;
    para_int(it,ms_down(it)) = para_sum(it,ms_down(it))* height_res;
    
    perp_int_ss(it,ms_down(it))= perp_sum_ss(it,ms_down(it))* height_res;
    para_int_ss(it,ms_down(it)) = para_sum_ss(it,ms_down(it))* height_res;
    
    int_MS(it,id_cb_lidar(it)) = (perp_sum(it,ms_down(it)) + ...
        para_sum(it,ms_down(it))) * height_res ;
    int_SS(it,id_cb_lidar(it)) = (perp_sum_ss(it,ms_down(it)) + ...
        para_sum_ss(it,ms_down(it))) * height_res ;
    
    for jz = 1:length(ms_down(it):ms_up(it))
    perp_int(it,ms_down(it)+jz) = perp_sum(it,ms_down(it)+jz)* height_res;
   
    para_int(it,ms_down(it)+jz) = para_sum(it,ms_down(it)+jz)* height_res;
    
    perp_int_ss(it,ms_down(it)+jz) = perp_sum_ss(it,ms_down(it)+jz)* height_res;
   
    para_int_ss(it,ms_down(it)+jz) = para_sum_ss(it,ms_down(it)+jz)* height_res;
    
    int_MS(it,ms_down(it)+jz) = (perp_sum(it,ms_down(it)+jz) + ...
        para_sum(it,ms_down(it)+jz)) * height_res ;
    
    int_SS(it,ms_down(it)+jz) = (perp_sum_ss(it,ms_down(it)+jz) + ...
        para_sum_ss(it,ms_down(it)+jz)) * height_res ;
    
    end 
        
end


% AS = ones(size(beta_atten_ss));
% 
%     for it = 1:nt
%     perp_int(it,:) = cumsum(perp(it),2).*height_res;
%     para_int(it,:) = cumsum(para(it),2).*height_res;
%     int_MS(it,:) = (perp_int(it,:) + para_int(it,:));
%     
%     perp_int_ss(it,:) = cumsum(perp_ss(it),2).*height_res;
%     para_int_ss(it,:) = cumsum(para_ss(it),2).*height_res;
%     int_SS(it,:) = (perp_int_ss(it,:) + para_int_ss(it,:));
%      
%     
%     end 
%  

AS = int_SS ./ int_MS;
AS(isnan(AS)) = 1 ;
 aa = ones(250,1);
 
AS = [AS aa]; 

%         

figure
loglog(beta_atten(100,1:end-1), height, '.-r', beta_atten_ss(100,2:end), height, '.-k')

save('correction_factor.mat', 'AS')
