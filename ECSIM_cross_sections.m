%% ECSIM Lidar inversion - plots only
% close all
max_z = max(z(id_ct_radar));
nt = size(Ze(:,:),1);
along_track = linspace(1,nt,nt);
%% Plot time series of measurements
TitleFigure=['ECSIM_time_series'];
figure('name', TitleFigure, 'NumberTitle','off', ...
    'units','centimeters','Position',[2 50 15 20]);
subplot(2,1,1)
Ze2= Ze;
Ze2(Ze2 ==0)=NaN;
p11 = pcolor(along_track,z,Ze2');
set(p11, 'EdgeColor', 'none');
set(gca,'Ydir','normal');
set(gca,'ylim',[0 max_z])
ylabel('z [m]')
hold on
l2 = plot(along_track,z(norm_down),'k','LineWidth',1.5);
hold on
l3 = plot(along_track,z(id_cb_lidar),'m','LineWidth',1.5);
set(gca, 'FontSize',12)
c=colorbar('southoutside');
c.Label.String = 'Z [mm^6/m^3]';
c.Label.FontSize = 12;
t = title('Radar Reflectivity Factor', 'FontSize',12,'FontWeight','normal');
set(t, 'horizontalAlignment', 'left')
set(t, 'units', 'normalized')
h1 = get(t, 'position');
set(t, 'position', [0 h1(2) h1(3)])
grid on
h_legend=legend([l2,l3], 'Start of the normalisation interval', 'Cloud base estimate');
set(h_legend,'FontSize',12);

subplot(2,1,2)
p12 = pcolor(along_track,z,Ba_MS');
set(p12, 'EdgeColor', 'none');
set(gca,'Ydir','normal');
set(gca,'ylim',[0 max_z])
hold on
plot(along_track,z(norm_down),'k','LineWidth',1.5)
hold on
plot(along_track,z(id_cb_lidar),'m','LineWidth',1.5)
ylabel('z [m]')
set(gca, 'FontSize',12)
t = title('Attenuated Backscatter Coefficient', 'FontSize',12,'FontWeight','normal');
set(t, 'horizontalAlignment', 'left')
set(t, 'units', 'normalized')
h1 = get(t, 'position');
set(t, 'position', [0 h1(2) h1(3)])
legend([l2,l3], 'Start of the normalisation interval', 'Cloud base estimate')
caxis([0 4*1.e-4])
c=colorbar('southoutside');
c.Label.String = 'ATB [m^{-1} sr^{-1}]';
% saveas(gcf,sprintf(TitleFigure),'epsc')
% export_fig(sprintf(TitleFigure), '-eps', '-transparent')



%% Plot time series of extinction
% TitleFigure=['Original and Retrieved Extinction'];
% figure('name', TitleFigure, 'NumberTitle','off', ...
%     'units','centimeters','Position',[2 50 15 20]);
% Ext2 = Ext;
% Ext2(Ext2==0)=NaN;
% subplot(2,1,1)
% p11 = pcolor(along_track,z,Ext2');
% set(p11, 'EdgeColor', 'none');
% set(gca,'Ydir','normal');
% set(gca,'ylim',[0 max_z])
% hold on
% plot(along_track,z(norm_down),'k','LineWidth',1.5)
% hold on
% plot(along_track,z(id_cb_lidar),'m','LineWidth',1.5)
% ylabel('z [m]')
% % hline = plot([tt], get(gca,'ylim') );
% % hline.Color = 'r';
% % hline.LineWidth = 1.5;    
% line('XData', [tt tt], 'YData', [0 max_z], 'LineWidth', 1.5, ...
%     'LineStyle', '-.', 'Color', 'r')
% % set(gca,'YTick',500:1000:max_z)
% % set(gca, 'YTickLabel', [0.5 1.5 2.5]); 
% set(gca, 'FontSize',8)
% c=colorbar('southoutside');
% caxis([0 0.06])
% c.Label.String = '[m^{-1}]';
% c.Label.FontSize = 8;
% % xlabel('Time [UTC]')
% t = title('Modeled Extinction from ECSIM', 'FontSize',10,'FontWeight','normal');
% set(t, 'horizontalAlignment', 'left')
% set(t, 'units', 'normalized')
% h1 = get(t, 'position');
% set(t, 'position', [0 h1(2) h1(3)])
% grid on
% 
% 
% subplot(2,1,2)
% p12 = pcolor(along_track,z,alpha');
% set(p12, 'EdgeColor', 'none');
% set(gca,'Ydir','normal');
% set(gca,'ylim',[0 max_z])
% ylabel('z [m]')
% set(gca, 'FontSize',8)
% t = title('Extinction Retrieved', 'FontSize',10,'FontWeight','normal');
% set(t, 'horizontalAlignment', 'left')
% set(t, 'units', 'normalized')
% h1 = get(t, 'position');
% set(t, 'position', [0 h1(2) h1(3)])
% caxis([0 0.06])
% c=colorbar('southoutside');
% c.Label.String = '[m^{-1}]';


%% Plot time series of extinction - the difference of the retrieved and original
TitleFigure=['Error_tau'];
figure('name', TitleFigure, 'NumberTitle','off', ...
    'units','centimeters','Position',[2 50 15 20]);
subplot(2,1,1)
tau_org2 = tau;
tau_org2(tau_org2==0)=NaN;
p12 = pcolor(along_track,z,tau_org2');
set(p12, 'EdgeColor', 'none');
set(gca,'Ydir','normal');
set(gca,'ylim',[0 600])
hold on
l1 = plot(along_track,z(norm_down),'k','LineWidth',1.5);
hold on
l2 = plot(along_track,z(id_cb_lidar),'m','LineWidth',1.5);
ylabel('z [m]')
set(gca, 'FontSize',12)
t = title('Cloud Optical Thickness', 'FontSize',12,'FontWeight','normal');
set(t, 'horizontalAlignment', 'left')
set(t, 'units', 'normalized')
h1 = get(t, 'position');
set(t, 'position', [0 h1(2) h1(3)])
grid on
legend([l1 l2],'Start of the normalisation interval', 'Cloud base estimate')
caxis([0 15])
c=colorbar('southoutside');
c.Label.String = '\tau';

subplot(2,1,2)

p11 = pcolor(along_track,z,err_alpha_MS_corr_plot');
set(p11, 'EdgeColor', 'none');
set(gca,'Ydir','normal');
set(gca,'ylim',[0 600])
hold on
plot(along_track,z(norm_down),'k','LineWidth',1.5)
hold on
plot(along_track,z(id_cb_lidar),'m','LineWidth',1.5)
ylabel('z [m]')
set(gca, 'FontSize',12)
c=colorbar('southoutside');
caxis([0 30])
c.Label.FontSize = 12;
c.Label.String = '%';
t = title('Retrieval Percent Error - with MS and RES correction', 'FontSize',12,'FontWeight','normal');
set(t, 'horizontalAlignment', 'left')
set(t, 'units', 'normalized')
h1 = get(t, 'position');
set(t, 'position', [0 h1(2) h1(3)])
grid on
% saveas(gcf,sprintf(TitleFigure),'epsc')
% export_fig(sprintf(TitleFigure), '-eps', '-transparent') 
% subplot(3,1,3)
% 
% p11 = pcolor(along_track,z,err_alpha_MS_plot');
% set(p11, 'EdgeColor', 'none');
% set(gca,'Ydir','normal');
% set(gca,'ylim',[0 400])
% hold on
% plot(along_track,z(norm_down),'k','LineWidth',1.5)
% hold on
% plot(along_track,z(id_cb_lidar),'m','LineWidth',1.5)
% ylabel('z [m]')
% set(gca, 'FontSize',10)
% c=colorbar('southoutside');
% caxis([0 30])
% c.Label.FontSize = 10;
% c.Label.String = '%';
% t = title('Retrieval Percent Error - with MS correction only', 'FontSize',12,'FontWeight','normal');
% set(t, 'horizontalAlignment', 'left')
% set(t, 'units', 'normalized')
% h1 = get(t, 'position');
% set(t, 'position', [0 h1(2) h1(3)])
% grid on

figure
p11 = pcolor(along_track,z,Ze');
set(p11, 'EdgeColor', 'none');
set(gca,'Ydir','normal');
%set(gca,'ylim',[0 400])
hold on
plot(along_track,z(norm_down),'k','LineWidth',1.5)
hold on
plot(along_track,z(id_cb_lidar),'m','LineWidth',1.5)
ylabel('z [m]')
set(gca, 'FontSize',10)
c=colorbar('southoutside');
%caxis([0 30])
c.Label.FontSize = 10;
%c.Label.String = '%';
t = title('Retrieval Percent Error - with MS correction only', 'FontSize',12,'FontWeight','normal');
set(t, 'horizontalAlignment', 'left')
set(t, 'units', 'normalized')
h1 = get(t, 'position');
set(t, 'position', [0 h1(2) h1(3)])
grid on
