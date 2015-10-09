% Interpolate signal with a 15 m resolution to 1 m resolution
% clear all;
% close;
format long
figure('NumberTitle','off', ...
    'Units', 'centimeters','Position',[2 50 45 15])

interp_factor = 5; %interpolation factor - how many measurements to add in between points
n_ext = [0.1 1 10 100]; % use 4 numbers or change subplot grid

for ie = 1: numel(n_ext);

% rho_atm = linspace(2.5*1.0e19,2.45*1.0e19,600)';
rho_atm = 2.55*1.0e19; %number concentration in cm-3
% rho should be in cm-3

rho_w=1.0e+6   ;% g/cm^3
% Height in km1
height=linspace(0.015,0.300,20)';
height_interp = linspace(0.015,0.300,80)';% all heights in km
% height_interp = linspace(0.015,0.300,80)';
%% Definie necessary parameters
cb = length(height)/2 ; %cloud base in the middle of the scene
% cb_coarse = length(height_coarse)/2 ;
cb_interp = length(height_interp)/2 ;

%% Calculate the molecular extinction and backscatter coefficients
% Calling ray_sigma_beta script, which provides Rayleigh scattering parameters
ray_sigma_beta
% Calculating molecular extinction and backscatter fields
alpha_m = rho_atm .* alpha_ray ;  % now in sr/cm
beta_m = rho_atm .* beta_ray   ; % now in sr/cm
height_res = height(2) - height(1);


% 1. Cloud droplet effective radius;
ext = zeros(length(height),1);
ext(1:cb-1) = alpha_m;
ext(cb:end) = alpha_m + n_ext(ie);
% ext = cumsum(ext);


% Calculate molecular optical thickness
tau_mol =  height_res .* cumsum(ext);
tau_mol(cb:end) = tau_mol(cb-1);
tau_cld = height_res .* cumsum(ext);

S = ones(length(height),1);
S(1:cb-1)=50;
S(cb:end)=20.0 ; 
S_interp = ones(length(height_interp),1);
S_interp(1:cb_interp-1)=50;
S_interp(cb_interp:end) = 20.0 ;

beta = ((ext./S) + beta_m) .* exp(-2.0 .* (tau_cld+tau_mol));
% beta = ext./S  .* exp(-2.0 .* tau_cld) ;
sd_beta_atten = 0.1 * beta ; 
P_sig = beta ./ height.^2;
P_sig_log = log(P_sig);
sd_P_sig = 0.2 .* P_sig;

% P_sig_interp_log = interp(P_sig_log,interp_factor);

P_sig_interp_log = interp1(height,P_sig_log,height_interp,'linear');
% 
P_sig_interp = exp(P_sig_interp_log);


height_res_int = height_interp(2)-height_interp(1);

% Find the peak of the attenuated backscatter
 Peak_P = max(P_sig(cb:end)) ;
 I_max  = find(P_sig >= Peak_P, 1, 'last' );
 [Peak_P_interp, I_max_interp] = max(P_sig_interp);
 

 
%% Calculate normalization value for each time step
% Define up and down limit for normalization
down = 5;
up = 7;

      
 %% Normalization height

 height_norm_1 = height(I_max+down) ;
 height_norm_2 = height(I_max+up) ; 
  
%  tau_mol = alpha_m * height_res;
 tau_m = 0;   
 int_S_beta = 0;

%% Calculating Pp, Pp_u, Pp_l and associated standard deviations in a loop

for j = 1: numel(height)
    
    tau_m = tau_m + alpha_m*height_res;
    int_S_beta = int_S_beta + (beta_m .* S(j)).*height_res;
    Pp(j) = S(j) .* P_sig(j) .* exp(-2.0.*int_S_beta) .* exp(2.*tau_m);
    sd_Pp(j) = S(j) .* sd_P_sig(j) .* exp(-2.0 .*int_S_beta) .* exp(2.*tau_m);   
          
end
clear int_S_beta tau_m
tau_m = 0;   
int_S_beta = 0;

for j = 1: numel(height_interp)
    
    int_S_beta = int_S_beta + (beta_m .* S_interp(j)) .* height_res_int;
    tau_m = tau_m + alpha_m*height_res_int;
    Pp_interp(j) = S_interp(j) .* P_sig_interp(j) .* exp(-2.0.*int_S_beta) .* exp(2.*tau_m);
%     sd_Pp_coarse(j) = S_coarse(j) .* sd_P_sig_(j) .* exp(-2.0 .*int_S_beta) .* exp(2.*tau_m);                            
        
end

height_o= (height_norm_1 + height_norm_2)/2;
io  = find(height <= height_o, 1, 'last' );
io_int = find(height_interp <= height_o, 1, 'last' );
  
alpha_o = ext(io) + S(io) .* beta_m;

x = log(beta(I_max+down:I_max+up));
y = height(I_max+down:I_max+up);



lm = fitlm(x,y);
% alpha_slope = -0.5 * (lm.Coefficients.Estimate(2)) + 1./x
alpha_slope = -0.5.* (diff(beta(I_max+down:I_max+up))./height_res) .* 1./beta((I_max+down:I_max+up-1));
alpha_slope2 = mean(-0.5.* (diff(x)./height_res));
% alpha_slope2 = alpha_slope2(round(length(beta_slope{ie})/2));
% alpha_o_int = ext(io) + S(io) .* beta_m;%
% alpha_o_coarse = ext_coarse(io_coarse);%
% alpha_o_int = ext_interp(io_int);%

%% Klett inversion  
RCS = Pp' .* height.^2 ;
RCS_interp = Pp_interp' .* height_interp.^2;

RCS_o = RCS(io) ;
% RCS_o_int = RCS_interp(io_int);

int_sig = zeros(1,numel(height));
int_sig_int = zeros(1,numel(height_interp));
int_sig(1) =  0.5 * RCS(1)*height_res;
int_sig_int(1) = 0.5 * RCS_interp(1)*height_res;

%   int_term[0]=0.5*P[0]*z[0]^2
%   for i=1,n-1 do begin
%      int_term[i]=int_term[i-1]+0.5*(RCS[i-1]+RCS[i])*res
%   endfor
% ;
%   int_term=(int_term-int_term[io])/RCS_o
 
for i=2:numel(height)-1;    
    int_sig(i) = int_sig(i-1) + 0.5*((RCS(i-1)+RCS(i))*height_res);    
end

for ic = 2:(numel(height_interp));
    int_sig_int(ic) = int_sig_int(ic-1) + 0.5 * (RCS_interp(ic-1) *height_res_int...
        + RCS_interp(ic)*height_res_int)  ;
end

  
int_S = (int_sig-int_sig(io)) ./ RCS_o ;
alphap = (RCS./RCS_o) ./ (1.0./alpha_o - 2.0 .* int_S') ;


int_S_int = (int_sig_int-int_sig_int(io_int)) ./ RCS_o;
alphap_int= (RCS_interp./RCS_o) ./ ((1.0./alpha_o) ...
    - (2.0 .* int_S_int')) ;

% 
% 
alpha = alphap - S.*beta_m ;
alpha_int = alphap_int - S_interp.*beta_m;
ret_alpha{ie} = alpha;
real_ext{ie} = ext;

%% Plot

subplot(1,4,ie)
p1 = plot(ext, height, '.-b');
hold on
p2 = plot(alpha,height, '*g');
hold on
p3 = plot(alpha_int,height_interp, '.m');
title([ num2str(n_ext (ie)) '/km set as extinction']);
hold on 
% normalisation = patch( ...
%     [height_org(I_max+down) height_org(I_max+down) height_org(I_max+up) height_org(I_max+up)],'r');
% set(normalisation,'FaceColor','b', ...
%     'FaceAlpha', 0.1)
% xlim([-0.1 10])
refline([0 height(I_max+down)])
refline([0 height(io)])





%% Plot
% 
subplot(1,4,ie)
p1 = plot(P_sig, height, 'd-k');
hold on
p2 = plot(P_sig_interp,height_interp, '.-r');

title(['Signal ' num2str(n_ext (ie)) '/km set as Ext']);

% subplot(1,3,ie)
% p1 = plot(Pp, height, 'd-b');
% hold on
% p2 = plot(Pp_interp,height_interp, '.-y');
% 
% title(['Signal ' num2str(n_ext (ie)) '/km set as Ext']);
accuracy(ie) = (abs(alpha_slope2-alpha_o)./alpha_o) * 100;
ac_max(ie) = max(accuracy(ie));
ac_min(ie) = min(accuracy(ie));

% acc_all = alpha./ext

end
rleg = legend(  ' True extinction', ...
   ' Retrieved (15 m)', ...
   ' Retrieved (1 m)');
legend('boxoff')

%% Simple check of accuracy of alpha_o
% ac_max
% ac_min
% 
% (ext-alpha)
error = (abs(ret_alpha{4} - real_ext{4})./real_ext{4}) * 100









 