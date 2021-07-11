% change the resolution of the grid to 15m - corresponding to the
% resolution of UV lidar

% Load data from ECSIM simulation
load('fire_27_merged_with_mask_v9a.mat')
load('fire_27_merged_with_mask_v9a_ss.mat')
load('true_fire_27_merged_with_mask_v9a.mat')
% load('correction_factor.mat')
% correction_SS

%% Select a specific part of the simulation
% Simulation runs from 0 to 250
start = 301;
stop = 800;
% Select needed varbles in the selected period
Ba_MS = lid_total(start:stop,1:end-8); % Original unit is 1/m/sr
depol_MS = lin_depol(start:stop,1:end-8);
% Ba_MS = awgn(Ba_MS,30,'measured');

Ba_SS = lid_total_ss(start:stop,1:end-8);
depol_SS = lin_depol_ss(start:stop,1:end-8);

% resize matrixes
Ba_MS_uv = resizem(Ba_MS, [500, 64]);
depol_MS_uv = resizem(depol_MS, [500, 64]);
Ba_SS_uv = resizem(Ba_SS, [500, 64]);
depol_SS_uv = resizem(depol_SS, [500, 64]);

% resize height
height = height(start:stop,1:end-8) .* 1.0e3;
z_uv = resizem(height, [500, 64]);
z_uv = z_uv(1, :);

% Divide multiple signal into parallel and perpendicular
para_ms = Ba_MS_uv ./ (depol_MS_uv + 1) ;
perp_ms = Ba_MS_uv - para_ms;

Ze = Ze(start:stop,1:end-8);
Ext = Ext(start:stop,9:end);
P = P(start:stop,1:end-8);
T = T(start:stop,1:end-8);
LWP = LWP(start:stop,1:end-8);


% sd_beta_atten = 0.2 .* Ba_MS ;
% 
% nz = size(lid_total(start:stop,1:end-8),2);   
% nt = size(lid_total(start:stop,1:end-8),1);
% along_track = linspace(1,nt,nt);
% z = height(1,1:end-8) .* 1.0e3;
% res = z(2) - z(1);  %Height resolution in meters