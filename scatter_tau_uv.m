%% Create scatter plots - divided by tau
for it = 1:nt
    for ir = 1:length(retrieved_cloud_uv{it})
    data_alpha_uv{ir} = alpha_ret_ms2ss_uv(:, retrieved_cloud_uv{it}(ir)) ; 
    data_ext_uv{ir} = ext_ray_uv(:,retrieved_cloud_uv{it}(ir));
    data_tau_uv{ir} = tau_uv(:,retrieved_cloud_uv{it}(ir));
    data_height_uv{ir} = z_uv(retrieved_cloud_uv{it}(ir)) - z_uv(retrieved_cloud_uv{it}(1));
    
    data_alpha_uv{ir}(data_alpha_uv{ir} ==0)=NaN;
    data_ext_uv{ir}(data_ext_uv{ir} ==0)=NaN;
    data_tau_uv{ir}(data_tau_uv{ir} ==0)=NaN;
    
    data_alpha_uv{ir}(data_alpha_uv{ir} < 0)=NaN;
    data_ext_uv{ir}(data_ext_uv{ir} < 0)=NaN;
    data_tau_uv{ir}(data_tau_uv{ir} < 0)=NaN;
    end
end

all_alpha = vertcat(data_alpha_uv{1},data_alpha_uv{2},data_alpha_uv{3}, ...
    data_alpha_uv{4},data_alpha_uv{5},data_alpha_uv{6},data_alpha_uv{7});
all_ext = vertcat(data_ext_uv{1},data_ext_uv{2},data_ext_uv{3}, ...
    data_ext_uv{4},data_ext_uv{5},data_ext_uv{6},data_ext_uv{7});
all_tau = vertcat(data_tau_uv{1},data_tau_uv{2},data_tau_uv{3}, ...
    data_tau_uv{4},data_tau_uv{5},data_tau_uv{6},data_tau_uv{7});

%% Define subplot net
sx = 2 ; % number of rows
sy = 2; %number of columns
n_bins = sx*sy; % must be divided by 2
max_tau = 2.25;
min_tau = 0.25 ;
bin_size =  ((max_tau - min_tau)/n_bins);

for i = 1:n_bins
 %clear id_tau
 id_tau = find(all_tau >= min_tau + bin_size*(i-1) & ...
         all_tau <= min_tau + bin_size*i);
 tau_id{i} = id_tau;
 % Divide data into bins based on the LWP
    for j = 1:length(tau_id{i})
        
      corr_alpha(i,j) = all_alpha(tau_id{i}(j));
      corr_ext(i,j) = all_ext(tau_id{i}(j));
      corr_tau(i,j) = all_tau(tau_id{i}(j));
    end
     
 % Calculate correlation coefficient for each scatter plot
  alpha_all{i} = (corr_alpha(i,:));
  ext_all{i} = (corr_ext(i,:));
  
  tau_all{i} = (corr_tau(i,:));
  tau_all{i} = tau_all{i}(~isnan(tau_all{i}));
 
    alpha_all{i}(alpha_all{i} ==0)=NaN;
    ext_all{i}(ext_all{i} ==0)=NaN;
    tau_all{i}(tau_all{i} ==0)=NaN;

end
% 
% 
TitleFigure=['scatter_tau'];
s = 55;
figure('name', TitleFigure, 'NumberTitle','off', ...
    'Units','centimeters','Position',[10 30 30 18]);
    for i = 1:n_bins
     cx(i) = subplot(sx,sy,i) ; 
     scat_fig_alpha(i) = scatter((alpha_all{i}),(ext_all{i}), s, ...
         tau_all{i} , 'fill');
     hold on
     legend('off')
     set(gca, 'FontSize',12)
     xlim([0 0.035])
     ylim([0 0.035])
     line('XData', [0 0.035], 'YData', [0 0.035], 'LineWidth', 1, ...
    'LineStyle', '-', 'Color', 'r')
     caxis([min_tau max_tau])
    a=get(gca,'XLim');
    x=max(a)-(max(a)/2.25);
    b=get(gca,'YLim');
    y=max(b)-(max(b)-min(b))/1.25;
    text(x,y,[...
    '\newline','\itA_{\alpha} = ', ...
    num2str(nanmean(abs(alpha_all{i}./ext_all{i})*100),'%2.2f%%') ...
    '\newline','\itE_{\alpha} = ',...
    num2str(nanmean((abs(alpha_all{i}-ext_all{i})./ext_all{i})*100),'%2.2f%%'), ...
    '\newline','n = ',...
    num2str(sum(~isnan(tau_all{i})),'%2.f')], ...
    'FontSize',12)
     title([num2str((min_tau + (bin_size*(i-1))),'%1.2f') ' < \tau < ',  ... 
       num2str((min_tau + (bin_size*i)),'%1.2f')],'FontSize',12,'FontWeight','normal')
    end
H=labelEdgeSubPlots('\alpha modelled [m^{-1} sr]','\alpha retrieved [m^{-1} sr]');
h2=colorbar;
set(h2, 'Position', [.925 .11 .0181 .81])
ylabel(h2,'\tau')
% h2.Label.String = '\tau';
h2.Label.FontSize = 12;
h2.Label.Rotation = 0; 
export_fig(sprintf(TitleFigure), '-eps', '-transparent')
saveas(gcf,sprintf(TitleFigure),'epsc')
% 
% 
