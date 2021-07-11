%% Lidar inversion plots
close all;

tt=10; %45 gives really nice profile

%% Corrections to retrieval
TitleFigure='retrieved_profile vs optical thickness';
figure('NumberTitle','off', ...
    'Units', 'centimeters','Position',[2 20 15 15])
p1 = plot( ...
    Ext(tt,id_cb_lidar(tt):norm_up(tt)), ...
    tau_org(tt,id_cb_lidar(tt):norm_up(tt)),'.-k');
hold on
p2 = plot(...
    alpha(tt,id_cb_lidar(tt):io(tt)), ...
    tau_org(tt,id_cb_lidar(tt):io(tt)),'.-b',...
    alpha_SS(tt,id_cb_lidar(tt):io(tt)), ...
    tau_org(tt,id_cb_lidar(tt):io(tt)), '.-g', ...
    alpha_corr(tt,id_cb_lidar(tt):io(tt)), ...
    tau_org(tt,id_cb_lidar(tt):io(tt)), '-.r', ...
    alpha_SS_corr(tt,id_cb_lidar(tt):io(tt)), ...
    tau_org(tt,id_cb_lidar(tt):io(tt)), '-.m', ...
    'LineWidth',1.5,...
    'MarkerSize',10);
% ylim([0 10])
t = title('Retrieved Extinction Coefficient');
set(t, 'horizontalAlignment', 'left')
set(t, 'units', 'normalized')
h1 = get(t, 'position');
set(t, 'position', [0 h1(2) h1(3)])
grid on
hold on
% xlim([0 0.05])
ylim([0 max(tau_org(tt,id_cb_lidar(tt):norm_up(tt)))])
line('XData', [0 0.05], 'YData', [tau_org(tt,norm_down(tt)) tau_org(tt,norm_down(tt))], ...
    'LineWidth', 1.0, 'LineStyle', '-.', 'Color', 'r')
% refline([0 tau_org(tt,norm_down(tt))]);
l = legend('True extinction', 'Retrieved extinction', ...
    'Retrieved extinction with MS correction', ...
    'Retrieved extinction with RES correction', ...
    'Retrieved extinction with MS & RES correction', ...
    'Begining of the normalisation invterval');
set(gca,'FontSize',12)
xlabel('\alpha [m^{-1} sr]')
ylabel('\tau ')
% saveas(gcf,sprintf(TitleFigure),'epsc')
% export_fig(sprintf(TitleFigure), '-eps', '-transparent')

TitleFigure='retrieved_profile vs height';
figure('NumberTitle','off', ...
    'Units', 'centimeters','Position',[20 20 15 15])
p1 = plot( ...
    (Ext(tt,id_cb_lidar(tt):norm_up(tt))), ...
    height(id_cb_lidar(tt):norm_up(tt)),'.-k');
hold on
p2 = plot(...
    (alpha(tt,id_cb_lidar(tt):io(tt))), ...
    height(id_cb_lidar(tt):io(tt)),'.-b',...
    (alpha_SS(tt,id_cb_lidar(tt):io(tt))), ...
    height(id_cb_lidar(tt):io(tt)), '.-g', ...
    (alpha_corr(tt,id_cb_lidar(tt):io(tt))), ...
    height(id_cb_lidar(tt):io(tt)), '.-r', ...
    (alpha_SS_corr(tt,id_cb_lidar(tt):io(tt))), ...
    height(id_cb_lidar(tt):io(tt)), '.-m', ...
    'LineWidth',1.5,...
    'MarkerSize',10);
% set(gca, 'XScale', 'log')
ylim([0 400])
t = title('Retrieved Extinction Coefficient');
set(t, 'horizontalAlignment', 'left')
set(t, 'units', 'normalized')
h1 = get(t, 'position');
set(t, 'position', [0 h1(2) h1(3)])
grid on
hold on
% xlim([0 0.05])
% ylim([0 max(height(id_cb_lidar(tt):norm_up(tt)))])
line('XData', [0 0.05], 'YData', [height(norm_down(tt)) height(norm_down(tt))], ...
    'LineWidth', 1.0, 'LineStyle', '-.', 'Color', 'r')
% refline([0 tau_org(tt,norm_down(tt))]);
l = legend('True extinction', 'Retrieved extinction', ...
    'Retrieved extinction with MS correction', ...
       'Retrieved extinction with RES correction', ...
    'Retrieved extinction with MS & RES correction', ...
    'Begining of the normalisation invterval');
set(gca,'FontSize',12)
xlabel('\alpha [m^{-1} sr]')
ylabel('[m] ')
% saveas(gcf,sprintf(TitleFigure),'epsc')
% export_fig(sprintf(TitleFigure), '-eps', '-transparent')

%% Simple check of accuracy of alpha_o
for i = 1:nt
    accuracy(i) = alpha_slope_SS(i)./Ext(i,io(i));
    bias(i) = alpha_slope_SS (i) -(Ext(i,io(i)) + S(it, io(it)) .* beta_m(it, io(it)));
end
mean_bias = mean(bias)
max_accu_alpha_o = max(accuracy(~isinf(accuracy)))
min_accu_alpha_o = min(accuracy(~isinf(accuracy)))
mean_accu_alpha_o = mean(accuracy(~isinf(accuracy)))
% 
%% Calculate statistical data and make tables

for it = 1:nt
retrieved_cloud{it} = id_cb_lidar(it):norm_down(it);
height_cloud{it} = height(retrieved_cloud{it});
end


for it = 1:nt
    for iih = 1:nh
      accu_alpha_plot(it,iih) = (alpha_SS_corr(it,iih) ./ Ext(it,iih));
      err_alpha_MS_corr_plot(it,iih) = (abs(alpha_SS_corr(it,iih) ...
     - Ext(it,iih))./ Ext(it,iih))*100; 
      err_alpha_MS_plot(it,iih) = (abs(alpha_SS(it,iih) ...
     - Ext(it,iih))./ Ext(it,iih))*100; 
 
    end  
    
for ih = 1:length(retrieved_cloud{it})
 accu_alpha(it,ih) = (alpha(it,retrieved_cloud{it}(ih))...
     ./ Ext(it,retrieved_cloud{it}(ih))) *100;
 err_alpha(it,ih) = (abs(alpha(it,retrieved_cloud{it}(ih)) ...
     - Ext(it,retrieved_cloud{it}(ih)))./ Ext(it,retrieved_cloud{it}(ih)))*100; 
 accu_alpha_sscat(it,ih) = (alpha_SS(it,retrieved_cloud{it}(ih)) ...
     ./ Ext(it,retrieved_cloud{it}(ih)))*100;
 err_alpha_sscat(it,ih) =(abs(alpha_SS(it,retrieved_cloud{it}(ih)) ...
     - Ext(it,retrieved_cloud{it}(ih)))./ Ext(it,retrieved_cloud{it}(ih)))*100; 
 
accu_alpha_c(it,ih) = (alpha_corr(it,retrieved_cloud{it}(ih)) ...
     ./ Ext(it,retrieved_cloud{it}(ih)))*100;
 err_alpha_c(it,ih) = (abs(alpha_corr(it,retrieved_cloud{it}(ih)) ...
     - Ext(it,retrieved_cloud{it}(ih)))./ Ext(it,retrieved_cloud{it}(ih)))*100; 
 

 accu_alpha_sscat_c(it,ih) = (alpha_SS_corr(it,retrieved_cloud{it}(ih)) ...
     ./ Ext(it,retrieved_cloud{it}(ih)))*100;
 err_alpha_sscat_c(it,ih) = (abs(alpha_SS_corr(it,retrieved_cloud{it}(ih)) ...
     - Ext(it,retrieved_cloud{it}(ih)))./ Ext(it,retrieved_cloud{it}(ih)))*100; 
 
 single_improvement(it,ih) =  (err_alpha(it,ih) - err_alpha_sscat(it,ih)) ;   
 bias_alpha(:,ih) = (sum(alpha(:,retrieved_cloud{it}(ih))) ...
     - sum(Ext(:,retrieved_cloud{it}(ih)))) ./nt ;
 RMSE(:,ih) = sqrt(mean((Ext(it,retrieved_cloud{it}(ih)) ...
     - alpha(it,retrieved_cloud{it}(ih))).^2));
 RMSE_ss(:,ih) = sqrt(mean((Ext(it,retrieved_cloud{it}(ih)) ...
     - alpha_SS(it,retrieved_cloud{it}(ih))).^2));
end

accu_alpha_plot(it,1:retrieved_cloud{it}-1)=NaN;
accu_alpha_plot(it,max(retrieved_cloud{it}):end)=NaN;
err_alpha_MS_corr_plot(it,1:retrieved_cloud{it}-1)=NaN;
err_alpha_MS_corr_plot(it,max(retrieved_cloud{it}):end)=NaN;

err_alpha_MS_plot(it,1:retrieved_cloud{it}-1)=NaN;
err_alpha_MS_plot(it,max(retrieved_cloud{it}):end)=NaN;
   
end
accu_alpha_plot(accu_alpha_plot ==0) = NaN;
err_alpha_MS_corr_plot(err_alpha_MS_corr_plot==0) = NaN;
err_alpha_MS_plot(err_alpha_MS_plot==0) = NaN;

accu_alpha(accu_alpha ==0) = NaN;
err_alpha(err_alpha==0)= NaN;
accu_alpha_sscat(accu_alpha_sscat==0)= NaN;
err_alpha_sscat(err_alpha_sscat==0)= NaN;
accu_alpha_c(accu_alpha_c==0)= NaN;
err_alpha_c(err_alpha_c==0)= NaN;
accu_alpha_sscat_c(accu_alpha_sscat_c==0)= NaN;
err_alpha_sscat_c(err_alpha_sscat_c==0)= NaN;

%% Latex table
fprintf('m from cb & ACU   & ERR   & ACU SS  & ERR SS   & ACU RES   & ERR RES& ACU SS RES & ERR SS RES& number\n')
for ii=1:length(retrieved_cloud{1})
    fprintf('%8.0f & %3.2f%% &  %3.2f%% &  %3.2f%% &  %3.2f%% &&  %3.2f%% & %3.2f%%&  %3.2f%% & %3.2f%% &  %4.0f \\\\ \n', ...
        height(retrieved_cloud{1}(ii)) - height(retrieved_cloud{1}(1)) , ...
        nanmean(accu_alpha(:,ii)), ...
        nanmean(err_alpha(:,ii)), ...
        nanmean(accu_alpha_sscat(:,ii)),...
        nanmean(err_alpha_sscat(:,ii)), ...
        nanmean(accu_alpha_c(:,ii)),...
        nanmean(err_alpha_c(:,ii)), ...
        nanmean(accu_alpha_sscat_c(:,ii)),...
        nanmean(err_alpha_sscat_c(:,ii)), ...
        sum(~isnan(err_alpha_sscat_c(:,ii))))
%         mean(single_improvement(:,ii)), ...
%         bias_alpha(:,ii), ...
%         nanmean(RMSE(ii)), ...
%         nanmean(RMSE_ss(ii)))

end
% 
% %% Plot profile with horizontal errorbar based on the bias of the measurements
% % figure('NumberTitle','off', ...
% %     'Units', 'centimeters','Position',[2 80 15 15])
% % hb1 = herrorbar(alpha_s(tt,id_cb_lidar(tt):norm_down(tt)), ...
% %     tau_org(tt,id_cb_lidar(tt):norm_down(tt)), ...
% %     err_alpha_sscat(tt,1:length(alpha_s(tt,id_cb_lidar(tt):norm_down(tt)))));
% % set(hb1, 'LineWidth',2)
% % %     bias_alpha(1:length(alpha_s(tt,id_cb_lidar(tt):norm_down(tt)))))
% % hold on
% % p12 = plot(Ext(tt,id_cb_lidar(tt):norm_down(tt)), ...
% %     tau_org(tt,id_cb_lidar(tt):norm_down(tt)), '-k','MarkerSize',8);
% % set(p12, 'LineWidth',2)
% % ylim([0 80])
% % xlabel('\alpha [m^{-1} sr]')
% % ylabel('\tau ')
% % t2 = title('Retrieved Extinction coefficient +/- bias');
% % set(t2, 'horizontalAlignment', 'left')
% % set(t2, 'units', 'normalized')
% % h1 = get(t2, 'position');
% % set(t2, 'position', [0 h1(2) h1(3)])
% % grid on
% 
% % figure('NumberTitle','off', ...
% %     'Units', 'centimeters','Position',[2 120 15 15])
% % s = 20;
% % sp1 = scatter(alpha_s(tt,id_cb_lidar(tt):norm_down(tt)),...
% %     Ext(tt,id_cb_lidar(tt):norm_down(tt))',s,tau_org(tt,id_cb_lidar(tt):norm_down(tt)),'filled');
% % % sp1 = plotmatrix(alpha_s(:,id_cb_lidar:norm_down),...
% % %     Ext(:,id_cb_lidar:norm_down));
% % ylim([0 0.05])
% % xlim([0 0.05])
% % xlabel('\alpha retrieved [m^{-1} sr]')
% % ylabel('\alpha modelled [m^{-1} sr] ')
% % t3 = title('Retrieved extinction vs Modelled extinction');
% % refline(1,0) 
% % lsline
% % colorbar
% % set(t3, 'units', 'normalized')
% 
% % Create scatter plots
% for it = 1:nt
%     for ir = 1:length(retrieved_cloud{it})
%     data_alpha{ir} = alpha_SS_corr( :, retrieved_cloud{it}(ir)) ; 
%     data_ext{ir} = Ext(:,retrieved_cloud{it}(ir));
%     data_tau{ir} = tau_org(:,retrieved_cloud{it}(ir));
%     data_height{ir} = height(retrieved_cloud{it}(ir)) - height(retrieved_cloud{it}(1));
%     
%     data_alpha{ir}(data_alpha{ir} ==0)=NaN;
%     data_ext{ir}(data_ext{ir} ==0)=NaN;
%     end
% end

% Plot scatter plot of retrieved and modelled extinction coefficient
% TitleFigure='scatter_retrieval';
% figure('name', TitleFigure, 'NumberTitle','off', ...
%     'Units','centimeters','Position',[10 30 25 15]);
% s = 25;
% for is = 1:6
%  cx(is) = subplot(2,3,is) ; 
%  scat_fig(is) = scatter(data_alpha{is},data_ext{is}, s,  data_tau{is} , 'fill');
%  hold on
% %  hline = refline(1,0);
% %  hline.Color = 'r';
% %  hline.LineWidth = 1 ;
%  caxis([min(min([data_tau{:}])) max(max([data_tau{:}]))])
%  xlim([0 0.035])
%  ylim([0 0.035])
%  line('XData', [0 0.04], 'YData', [0 0.04], 'LineWidth', 1, ...
%     'LineStyle', '-', 'Color', 'r')
%  a=get(gca,'XLim');
%  x=max(a)-(max(a)/1.25);
%  b=get(gca,'YLim');
%  y=max(b)-(max(b)-min(b))/12;
%  text(x,y,[...%'\itr = ',num2str(correlations_reff_beta{i}(2),'%.2f'),...
%     ...%'\newline','\itr^2 = ',num2str(mdl_lidar_reff{i}.Rsquared.Ordinary,'%.2f'),...
%     '\newline','\itA_{\alpha} = ',num2str(nanmean(accu_alpha_sscat_c(:,is)),'%3.2f%%') ...
%     '\newline','\itE_{\alpha} = ',num2str(nanmean(err_alpha_sscat_c(:,is)),'%3.2f%%')], ...
%     'FontSize',10);
% title([num2str(data_height{is}) ' m from cloud base'],'FontSize',10,'FontWeight','normal')
% end
% H=labelEdgeSubPlots('\alpha modelled [m^{-1} sr]','\alpha retrieved [m^{-1} sr]');
% h2=colorbar;
% set(h2, 'Position', [.925 .11 .0181 .81])
% ylabel(h2,'\tau')
% % h2.Label.String = '\tau';
% h2.Label.FontSize = 10;
% h2.Label.Rotation = 0; 
% export_fig(sprintf(TitleFigure), '-eps', '-transparent')
% % saveas(gcf,sprintf(TitleFigure),'epsc')
% 
% % Create scatter plots
% for it = 1:nt
%     for ir = 1:length(retrieved_cloud{it})
%     data_alpha{ir} = alpha_SS_corr( :, retrieved_cloud{it}(ir)) ; 
%     data_ext{ir} = Ext(:,retrieved_cloud{it}(ir));
%     data_tau{ir} = tau_org(:,retrieved_cloud{it}(ir));
%     data_height{ir} = height(retrieved_cloud{it}(ir)) - height(retrieved_cloud{it}(1));
%     
%     data_alpha{ir}(data_alpha{ir} ==0)=NaN;
%     data_ext{ir}(data_ext{ir} ==0)=NaN;
%     end
% end
% 
% % Plot scatter plot of retrieved and modelled extinction coefficient
% TitleFigure='scatter_retrieval';
% figure('name', TitleFigure, 'NumberTitle','off', ...
%     'Units','centimeters','Position',[10 30 25 15]);
% s = 25;
% for is = 1:6
%  cx(is) = subplot(2,3,is) ; 
%  scat_fig(is) = scatter(data_alpha{is},data_ext{is}, s,  data_tau{is} , 'fill');
%  hold on
% %  hline = refline(1,0);
% %  hline.Color = 'r';
% %  hline.LineWidth = 1 ;
%  caxis([min(min([data_tau{:}])) max(max([data_tau{:}]))])
%  xlim([0 0.035])
%  ylim([0 0.035])
%  line('XData', [0 0.04], 'YData', [0 0.04], 'LineWidth', 1, ...
%     'LineStyle', '-', 'Color', 'r')
%  a=get(gca,'XLim');
%  x=max(a)-(max(a)/1.25);
%  b=get(gca,'YLim');
%  y=max(b)-(max(b)-min(b))/12;
%  text(x,y,[...%'\itr = ',num2str(correlations_reff_beta{i}(2),'%.2f'),...
%     ...%'\newline','\itr^2 = ',num2str(mdl_lidar_reff{i}.Rsquared.Ordinary,'%.2f'),...
%     '\newline','\itA_{\alpha} = ',num2str(nanmean(accu_alpha_sscat_c(:,is)),'%3.2f%%') ...
%     '\newline','\itE_{\alpha} = ',num2str(nanmean(err_alpha_sscat_c(:,is)),'%3.2f%%')], ...
%     'FontSize',10);
% title([num2str(data_height{is}) ' m from cloud base'],'FontSize',10,'FontWeight','normal')
% end
% H=labelEdgeSubPlots('\alpha modelled [m^{-1} sr]','\alpha retrieved [m^{-1} sr]');
% h2=colorbar;
% set(h2, 'Position', [.925 .11 .0181 .81])
% ylabel(h2,'\tau')
% % h2.Label.String = '\tau';
% h2.Label.FontSize = 10;
% h2.Label.Rotation = 0; 
% export_fig(sprintf(TitleFigure), '-eps', '-transparent')

%% Signal divided by signal
% figure('NumberTitle','off', ...
%     'Units', 'centimeters','Position',[2 20 15 15])
% p2 = plot((P_sig(tt,id_cb_lidar(tt):norm_down(tt))./P_sig_SS(tt,id_cb_lidar(tt):norm_down(tt))), ...
%     tau_org(tt,id_cb_lidar(tt):norm_down(tt)),'-.r', ...
%     (P_sig(tt,id_cb_lidar(tt):norm_down(tt))./P_sig_corr(tt,id_cb_lidar(tt):norm_down(tt))), ...
%     tau_org(tt,id_cb_lidar(tt):norm_down(tt)),'-.k', ...
%     (P_sig_corr(tt,id_cb_lidar(tt):norm_down(tt))./P_sig_SS_corr(tt,id_cb_lidar(tt):norm_down(tt))), ...
%     tau_org(tt,id_cb_lidar(tt):norm_down(tt)),'-.g', ...
%     (P_sig(tt,id_cb_lidar(tt):norm_down(tt))./P_sig_SS_corr(tt,id_cb_lidar(tt):norm_down(tt))), ...
%     tau_org(tt,id_cb_lidar(tt):norm_down(tt)),'-.b', ...
%     'LineWidth',1.5, ...
%     'MarkerSize',10);
% % ylim([0 10])
% t = title('P_{sig}/P_{sig_{single}}');
% set(t, 'horizontalAlignment', 'left')
% set(t, 'units', 'normalized')
% h1 = get(t, 'position');
% set(t, 'position', [0 h1(2) h1(3)])
% grid on
% legend('P / P_{SS}', 'P / P_{att}', 'P_{att} / P_{SS} _{att}' , ...
%     'P / P_{SS} _{att}')
% xlabel(' ')
% ylabel('\tau ')

%% Depolarisation ratio
% figure('NumberTitle','off', ...
%     'Units', 'centimeters','Position',[2 20 15 15])
% p3 = plot(depol_ratio(tt,id_cb_lidar(tt):norm_down(tt)), ...
%     tau_org(tt,id_cb_lidar(tt):norm_down(tt)),'-.r','LineWidth',1.5, ...
%     'MarkerSize',10);
% % ylim([0 10])
% t = title('Depolarisation ration');
% set(t, 'horizontalAlignment', 'left')
% set(t, 'units', 'normalized')
% h1 = get(t, 'position');
% set(t, 'position', [0 h1(2) h1(3)])
% grid on
% % legend('Retrieved extinction', '"True" extinction', 'Retrieved extinction with multiple scattering correction')
% xlabel('depolarisation ratio')
% ylabel('\tau ')

%% Plot alpha slope
% t1 = 1;
% TitleFigure='Slope_method';
% figure('name', TitleFigure, 'NumberTitle','off', ...
%     'Units','centimeters','Position',[10 30 15 15]);
% p1 = plot( ...  
%     Ext(t1,1:end-1), height(1:end-1),'.-r',...
%     alpha_o_p(t1,:), height(1:end-1), '.-k', ...
%     'LineWidth',1.5,...
%     'MarkerSize',10);
% % ylim([0 10])
% t = title('Slope Method Extinction Retrieval');
% set(t, 'horizontalAlignment', 'left')
% set(t, 'units', 'normalized')
% h1 = get(t, 'position');
% set(t, 'position', [0 h1(2) h1(3)])
% grid on
% hold on
% % xlim([0-mean(Ext(t1,1:end-1)) 0.06])
% ylim([0 height(norm_up(t1) +5)])
% line('XData', [0 max(Ext(t1,1:end-1))], 'YData', [height(io(t1)) height(io(t1))], 'LineWidth', 1.5, ...
%      'Color', 'b')
% line('XData', [0 max(Ext(t1,1:end-1))], 'YData', [height(id_cb_lidar(t1)) height(id_cb_lidar(t1))], 'LineWidth', 1.5, ...
%      'Color', 'g')
% legend('"True" extinction', 'Slope method retrieval', 'Normalisation height', 'Cloud base height')
% xlabel('\alpha [m^{-1} sr]')
% ylabel('Height [m] ')
% set(gca,'FontSize',12)
% % saveas(gcf,sprintf(TitleFigure),'epsc')
% export_fig(sprintf(TitleFigure), '-eps', '-transparent')
% 
% TitleFigure='Signal';
% figure('name', TitleFigure, 'NumberTitle','off', ...
%     'Units','centimeters','Position',[10 30 15 15]);
% p1 = plot( ...  
%     beta_atten_SS(t1,1:end-1), height(1:end-1),'.-r',...
%     'LineWidth',1.5,...
%     'MarkerSize',10);
% % ylim([0 10])
% t = title('Attenuated Backscatter Coefficient Profile');
% set(t, 'horizontalAlignment', 'left')
% set(t, 'units', 'normalized')
% h1 = get(t, 'position');
% set(t, 'position', [0 h1(2) h1(3)])
% grid on
% hold on
% % xlim([0 max(beta_atten_SS(t1,1:end-1))+mean(beta_atten_SS(t1,1:end-1))])
% % ylim([0 height(norm_up(t1) +50)])
% % line('XData', [0 max(beta_atten_SS(t1,1:end-1))], 'YData', [height(norm_up(t1)) height(norm_up(t1))], 'LineWidth', 1.5, ...
% %      'Color', 'b')
% % line('XData', [0 max(beta_atten_SS(t1,1:end-1))], 'YData', [height(norm_down(t1)) height(norm_down(t1))], 'LineWidth', 1.5, ...
% %     'Color', 'b')
% legend('ATB', 'Normalisation interval boundaries')
% xlabel('ATB [m^{-1} sr^{-1}]')
% ylabel('Height [m] ')
% set(gca,'FontSize',12)
% % saveas(gcf,sprintf(TitleFigure),'epsc')
% export_fig(sprintf(TitleFigure), '-eps', '-transparent')
