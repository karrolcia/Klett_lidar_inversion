format long
%% Version 2.0
% Correction for the multiple scattering is here calculated from
% B_ss = M * B_ms + I_ms * (dM/dz)
% where B_ss is single scattering attenuated backscatter,
% M is the correction factor calculated accordint to Hu et al., 2006
% I_ms is the integrated multiple scattering attenuated backscatter
% dM/dz is the derivative of the corrections factor

clearvars; close all;

%% Change the resolution of the grid to 15m - corresponding to the
% resolution of UV lidar

% Load data from ECSIM simulation
load('fire_27_merged_with_mask_v9a.mat')
load('fire_27_merged_with_mask_v9a_ss.mat')
load('true_fire_27_merged_with_mask_v9a.mat')
% load('correction_factor.mat')
% correction_SS

%% Select a specific part of the simulation
% Simulation runs from 0 to 250
start = 350;
stop = 800;
% Select needed varbles in the selected period
Ba_MS = lid_total(start:stop,1:end-8); % Original unit is 1/m/sr
depol_MS = lin_depol(start:stop,1:end-8);

Ba_SS = lid_total_ss(start:stop,1:end-8);
depol_SS = lin_depol_ss(start:stop,1:end-8);

height = height(start:stop,1:end-8) .* 1.0e3;
z = height(1,:);


%% Extract other needed information
Ze = Ze(start:stop,1:end-8);
Ze_uv = resizem(Ze, [500, 64]);
Ext = Ext(start:stop,9:end);
P = P(start:stop,1:end-8);
T = T(start:stop,1:end-8);
LWP = LWP(start:stop,1:end-8);

%% Resize matrixes
Ba_MS_uv = resizem(Ba_MS, [length(start:stop), 64]);
depol_MS_uv = resizem(depol_MS, [length(start:stop), 64]);
Ba_SS_uv = resizem(Ba_SS, [length(start:stop), 64]);
z_uv = resizem(height, [length(start:stop), 64]);
z_uv = z_uv(1, :);
Ze_uv = resizem(Ze, [length(start:stop), 64]);
Ext_uv = resizem(Ext, [length(start:stop), 64]);
P_uv = resizem(P, [length(start:stop), 64]);
T_uv = resizem(T, [length(start:stop), 64]);

% Set normalisation height (number of bins over cloud base)
nh = 40;
nh_uv = 6;

[alpha_ret_ss, ...
          alpha_ret_ss_nc, ...
          alpha_ret_ms, ...
          alpha_ret_ms_nc, ...
          alpha_ret_ms2ss_nc, ...
          alpha_ret_ms2ss, ...
          id_cb_lidar, ...
          norm_down, ...
          norm_up, ...
          id_ct_radar, tau, ext_ray, io, alpha_ss_slope, ...
          alpha_ss_slope_corr] = lidar_inversion_function(Ba_MS, ...
                                              depol_MS, ...
                                              Ba_SS, ...
                                              z, ...
                                              Ze, ...
                                              Ext, ...
                                              P, ...
                                              T, nh);
                                          
                     


[alpha_ret_ss_uv, ...
          alpha_ret_ss_nc_uv, ...
          alpha_ret_ms_uv, ...
          alpha_ret_ms_nc_uv, ...
          alpha_ret_ms2ss_nc_uv, ...
          alpha_ret_ms2ss_uv, ...
          id_cb_lidar_uv, ...
          norm_down_uv, ...
          norm_up_uv, ...
          id_ct_radar_uv, tau_uv, ext_ray_uv, io_uv, alpha_ss_slope_uv, ...
          alpha_ss_slope_corr_uv] = lidar_inversion_function(Ba_MS_uv, ...
                                              depol_MS_uv, ...
                                              Ba_SS_uv, ...
                                              z_uv, ...
                                              Ze_uv, ...
                                              Ext_uv, ...
                                              P_uv, ...
                                              T_uv, nh_uv);                                          
                                          
                                          
                                          


%
lida_inversion_tables_15
%lida_inversion_tables_2_5

ECSIM_cross_sections_uv
scatter_tau_uv 
plots_uv


%% Error calculation
%% Simple check of accuracy of alpha_o
nt = size(Ba_MS(:,:),1);
for i = 1:nt
    accuracy_ss(i) = alpha_ss_slope_uv(i,io_uv(i))./Ext_uv(i,io_uv(i));
    accuracy_ms(i) = alpha_ss_slope_corr_uv(i,io_uv(i))./Ext_uv(i,io_uv(i));
    bias(i) = alpha_ss_slope_uv(i) - (Ext_uv(i,io_uv(i)));
end
mean_bias = mean(bias)
max_accu_alpha_o = max(accuracy_ss(~isinf(accuracy_ss)))
min_accu_alpha_o = min(accuracy_ss(~isinf(accuracy_ss)))
mean_accu_alpha_o = mean(accuracy_ss(~isinf(accuracy_ss)))

max_accu_alpha_o_ms = max(accuracy_ms(~isinf(accuracy_ms)))
min_accu_alpha_o_ms = min(accuracy_ms(~isinf(accuracy_ms)))
mean_accu_alpha_o_ms = mean(accuracy_ms(~isinf(accuracy_ms)))
% 