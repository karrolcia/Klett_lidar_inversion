% Interpolate signal with a 15 m resolution to 1 m resolution
% clear all;
% close;
figure('NumberTitle','off', ...
    'Units', 'centimeters','Position',[2 20 15 15])
load('beta.mat')
% att_beta in sr/m
% height in m



rho_atm = 2.55*1.0e19; %number concentration in cm-3
% rho should be in cm-3

rho_w=1.0e+6   ;% g/cm^3
% Height in km1

% height_interp = linspace(15,300,300)'./1000;% all heights in km
%% Definie necessary parameters

%% Calculate the molecular extinction and backscatter coefficients
% Calling ray_sigma_beta script, which provides Rayleigh scattering parameters
ray_sigma_beta
% Calculating molecular extinction and backscatter fields
alpha_m = rho_atm .* alpha_ray * 1.0e-3;  % now in sr/m
beta_m = rho_atm .* beta_ray * 1.0e-3; % now in sr/m
height_res = height(2) - height(1);


P_sig = att_beta  ./ height.^2;
% P_sig_log = log(P_sig);
sd_P_sig = 0.2 .* P_sig;


Peak_P = max(P_sig(1:end)) ;
I_max  = find(P_sig >= Peak_P, 1, 'last' );
 
S = ones(length(height),1);
S(1:I_max-1)=50;
S(I_max:end)=18.0 ; 

%% Calculate normalization value for each time step
% Define up and down limit for normalization
down = 6;
up = 9;

lower_bound = I_max+down;
upper_bound = I_max+up;
      
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


height_o= (height_norm_1 + height_norm_2)/2;
io  = find(height <= height_o, 1, 'last' );
alpha_o = ext(io) + S(io) .* beta_m;

% Calculate slope of the signal for the inversion

% Calculating slope of the signal  
x = log(att_beta(I_max+down:I_max+up));
y = height(I_max+down:I_max+up);

  
 [poly_coeff , struct_S, mu] = polyfit(x,y,1);
 
  slope = poly_coeff(1);

  [yfit, d_alpha_o] = polyval(poly_coeff,x,struct_S);
  yresid = y - yfit;
  SSresid = sum(yresid.^2);
  SStotal = (length(y)-1) * var(y);
  rsq = 1 - SSresid./SStotal;

  
  lm = fitlm(x,y);
  alpha_o_p = -0.5 * (diff(log(att_beta(I_max+down:I_max+up)))./height_res) +  beta_m * S(io);
alpha_slope = alpha_o_p(2); % Take the second value as it corresponds to io
%% Klett inversion  
RCS = Pp .* height.^2 ;

RCS_o = RCS(io) ;
% RCS_o_int = RCS_interp(io_int);

int_sig = zeros(1,numel(height));
int_sig(1) = 0.5 * (RCS(1) * height_res);
 
for i=2:numel(height)-1;    
    int_sig(i) = int_sig(i-1) + (0.5 * (RCS(i-1)+RCS(i))) *height_res;    
end
  
int_S = (int_sig-int_sig(io)) / RCS_o ;
alphap = (RCS./RCS_o) ./ (1.0./alpha_slope - 2.0 .* int_S) ;


% 
% 
alpha = alphap';

%% Plot
subplot(1,1,1)
p1 = plot(ext(1:30), height(1:30), '.-b');
hold on
p2 = plot(alpha(1:30),height(1:30), '*-g');
hold on
% normalisation = patch( ...
%     [height_org(I_max+down) height_org(I_max+down) height_org(I_max+up) height_org(I_max+up)],'r');
% set(normalisation,'FaceColor','b', ...
%     'FaceAlpha', 0.1)
% xlim([-0.1 5])
refline([0 height(I_max+down)])
refline([0 height(io)])






%% Plot
% 
p3 = plot(P_sig(1:30), height(1:30), 'd-k');
title('Signal from ECSIM');

% subplot(1,3,ie)
% p1 = plot(Pp, height, 'd-b');
% hold on
% p2 = plot(Pp_interp,height_interp, '.-y');
% 
% title(['Signal ' num2str(n_ext (ie)) '/km set as Ext']);

rleg2 = legend(  ' True extinction', ...
   ' Retrieved (15 m)');
legend('boxoff')


figure('NumberTitle','off', ...
    'Units', 'centimeters','Position',[20 20 15 15])

her = herrorbar(alpha(1:30), height(1:30), ext(1:30)'-alpha(1:30));

hold on
plot(ext(1:30), height(1:30));
hold on
ref1 = refline([0 height(I_max+down)]);
ref2 = refline([0 height(I_max+up)]);
set(ref1, 'Color', 'black')
set(ref2, 'Color', 'black')
title('Retrieved extinction and Original extinction comparison');
rleg = legend(  'Original extinction-retrieved extinction', ...
    ' Retrieved extinction', 'Original extinction');
alpha_o
alpha_o_p
e = ext(1:30)'-alpha(1:30);
per_diff = e./ext(1:30)';
