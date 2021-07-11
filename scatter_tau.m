%% Create scatter plots - divided by tau
for it = 1:nt
    for ir = 1:length(retrieved_cloud{it})
    data_alpha{ir} = alpha_ret_ms2ss(:, retrieved_cloud{it}(ir)) ; 
    data_ext{ir} = ext_ray(:,retrieved_cloud{it}(ir));
    data_tau{ir} = tau(:,retrieved_cloud{it}(ir));
    data_height{ir} = z(retrieved_cloud{it}(ir)) - height(retrieved_cloud{it}(1));
    
    data_alpha{ir}(data_alpha{ir} ==0)=NaN;
    data_ext{ir}(data_ext{ir} ==0)=NaN;
    data_tau{ir}(data_tau{ir} ==0)=NaN;
    end
end

all_alpha = vertcat(data_alpha{1},data_alpha{2},data_alpha{3}, ...
    data_alpha{4},data_alpha{5},data_alpha{6},data_alpha{7});
all_ext = vertcat(data_ext{1},data_ext{2},data_ext{3}, ...
    data_ext{4},data_ext{5},data_ext{6},data_ext{7});
all_tau = vertcat(data_tau{1},data_tau{2},data_tau{3}, ...
    data_tau{4},data_tau{5},data_tau{6},data_tau{7});

%% Define subplot net
sx = 2; % number of rows
sy = 2; %number of columns
n_bins = sx*sy; % must be divided by 2
max_tau = 3;
min_tau = 0.3;
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


end
% 
% 
TitleFigure=['scatter_tau'];
s = 35;
figure('name', TitleFigure, 'NumberTitle','off', ...
    'Units','centimeters','Position',[10 30 25 15]);
    for i = 1:n_bins
     cx(i) = subplot(sx,sy,i) ; 
     scat_fig_alpha(i) = scatter(alpha_all{i},(ext_all{i}), s, ...
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
    num2str(nanmean((alpha_all{i}./ext_all{i})*100),'%2.2f%%') ...
    '\newline','\itE_{\alpha} = ',...
    num2str(nanmean((abs(alpha_all{i}-ext_all{i})./ext_all{i})*100),'%2.2f%%')], ...
    'FontSize',12);
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
% export_fig(sprintf(TitleFigure), '-eps', '-transparent')
% % saveas(gcf,sprintf(TitleFigure),'epsc')
% 
% 
% % % Plot scatter plot of retrieved and modelled extinction coefficient
% % TitleFigure='scatter_retrieval depending on tau';
% % figure('name', TitleFigure, 'NumberTitle','off', ...
% %     'Units','centimeters','Position',[10 30 25 15]);
% % s = 25;
% % for is = 1:size(
% %  cx(is) = subplot(2,3,is) ; 
% %  scat_fig(is) = scatter(data_alpha{is},data_ext{is}, s,  data_tau{is} , 'fill');
% %  hold on
% % %  hline = refline(1,0);
% % %  hline.Color = 'r';
% % %  hline.LineWidth = 1 ;
% %  caxis([min(min([data_tau{:}])) max(max([data_tau{:}]))])
% %  xlim([0 0.035])
% %  ylim([0 0.035])
% %  line('XData', [0 0.04], 'YData', [0 0.04], 'LineWidth', 1, ...
% %     'LineStyle', '-', 'Color', 'r')
% %  a=get(gca,'XLim');
% %  x=max(a)-(max(a)/1.25);
% %  b=get(gca,'YLim');
% %  y=max(b)-(max(b)-min(b))/12;
% %  text(x,y,[...%'\itr = ',num2str(correlations_reff_beta{i}(2),'%.2f'),...
% %     ...%'\newline','\itr^2 = ',num2str(mdl_lidar_reff{i}.Rsquared.Ordinary,'%.2f'),...
% %     '\newline','\itA_{\alpha} = ',num2str(nanmean(accu_alpha_sscat_c(:,is)),'%3.2f%%') ...
% %     '\newline','\itE_{\alpha} = ',num2str(nanmean(err_alpha_sscat_c(:,is)),'%3.2f%%')], ...
% %     'FontSize',10);
% % title([num2str(data_height{is}) ' m from cloud base'],'FontSize',10,'FontWeight','normal')
% % end
% % H=labelEdgeSubPlots('\alpha modelled [m^{-1} sr]','\alpha retrieved [m^{-1} sr]');
% % h2=colorbar;
% % set(h2, 'Position', [.925 .11 .0181 .81])
% % ylabel(h2,'\tau')
% % % h2.Label.String = '\tau';
% % h2.Label.FontSize = 10;
% % h2.Label.Rotation = 0; 
% % export_fig(sprintf(TitleFigure), '-eps', '-transparent')
% % % saveas(gcf,sprintf(TitleFigure),'epsc')
% 
