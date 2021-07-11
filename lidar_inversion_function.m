%% Lidar inversion function 

function [alpha_ret_ss, ...
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
                                              T, nh)


% 1. Divide signal into paralel and perpendicular compomnent based on the
% depolarisation
para_ms = Ba_MS ./ (depol_MS + 1) ;
perp_ms = Ba_MS - para_ms;

% Standard deviation of the Attenuated backscatter
sd_beta_atten = 0.2 .* Ba_MS ;

% Extract matrixes size for calculations

% the normalisation interval
nz = size(Ba_MS(:,:),2);   
nt = size(Ba_MS(:,:),1);
along_track = linspace(1,nt,nt);
res = z(2) - z(1);  %Height resolution in meters

%% Estimating cloud base height and defining normalisation interval
 %[Peak_P, I_max] = max(Ba_MS, [], 2) ;
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
     
    id_cb_lidar(it) = find (perp_ms(it,:) >= (Peak_perp_ms(it)/10), ...
                                                1,'first'); 
    id_ct_radar(it) = find((Ze(it,:)) >= 1.e-4, 1, 'last');

    norm_down(it) = id_cb_lidar(it) + nh; % Normalisation height 
    % 30*2.5 = 75 m above the cloud base
    norm_up(it) = norm_down(it) + 5;    % Higher limit of 
    % the normalisation height - normalize over 5 bins

    z_norm_1(it) = z(norm_down(it));
    z_norm_2(it) = z(norm_up(it)) ;

    z_o(it) = (z_norm_1(it) + z_norm_2(it))/2;
    io(it)  = find(z<= z_o(it), 1, 'last' );

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
      
    int_MS(it,ms_down(it):ms_up(it)) = ...
        (perp_int(it,ms_down(it):ms_up(it)) + ...
        para_int(it,ms_down(it):ms_up(it)));  
end
depol_ratio = perp_int ./ para_int ;
correction_factor = ((1 - depol_ratio) ./ (1 + depol_ratio)).^2;
correction_factor(isnan(correction_factor)) = 1 ;

dM = correction_factor;

for it= 1:nt
 for j = 2:length(z)-1
     dM(it,j) = dM(it,j+1) - dM(it,j);
 end
end
dMdz = dM ./ res;

Ba_SS_corr = Ba_SS_corr./res;
Ba_SS_corr2 = correction_factor .* Ba_MS + int_MS .* dMdz;

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
    FK = 1.225e-13*wavelen^4 - 3.911e-10*wavelen^3 + ...
        4.6100e-7*wavelen^2 - 2.410e-4*wavelen + 1.095 ;  
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
% Calculate tau_ray and int_beta_ray for SS
  tau=zeros(size(Ba_SS));
  tau_ray=zeros(size(Ba_SS));
  int_beta_ray=zeros(size(Ba_SS));
    
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

  
%% Calucalte single and multiple attenuated backscatter 
%  corrected for the Reayleigh   
  S=16.0;
  Ba_msp=S.*Ba_MS.*exp(2.*tau_ray).*exp(-2.*S.*int_beta_ray);
  Ba_ssp=S.*Ba_SS.*exp(2.*tau_ray).*exp(-2.*S.*int_beta_ray);
  Ba_ssp_corr = S.*Ba_SS_corr2.*exp(2.*tau_ray).*exp(-2.*S.*int_beta_ray);
  
%% Perform signal inversion
alpha_ret_ss = zeros(size(Ba_SS));
alpha_ret_ss_nc = zeros(size(Ba_SS));
alpha_ret_ms = zeros(nt,nz);
alpha_ret_ms_nc = zeros(nt,nz);
alpha_ret_ms2ss_nc = zeros(nt,nz);
alpha_ret_ms2ss = zeros(nt,nz);

 for it = 1:nt
     ext_ray(it,:) = extp(it,:)+S.*beta_ray(it,:);
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
 
alpha_ret_ss(alpha_ret_ss < 0) = NaN;
alpha_ret_ss_nc(alpha_ret_ss_nc < 0) = NaN;
alpha_ret_ms(alpha_ret_ms< 0) = NaN;
alpha_ret_ms_nc(alpha_ret_ms_nc< 0) = NaN;
alpha_ret_ms2ss_nc(alpha_ret_ms2ss_nc <0) = NaN;
alpha_ret_ms2ss(alpha_ret_ms2ss< 0) = NaN;

alpha_ret_ss(alpha_ret_ss > 2) = NaN;
alpha_ret_ss_nc(alpha_ret_ss_nc > 2) = NaN;
alpha_ret_ms(alpha_ret_ms > 2) = NaN;
alpha_ret_ms_nc(alpha_ret_ms_nc > 2) = NaN;
alpha_ret_ms2ss_nc(alpha_ret_ms2ss_nc > 2) = NaN;
alpha_ret_ms2ss(alpha_ret_ms2ss > 2) = NaN;
