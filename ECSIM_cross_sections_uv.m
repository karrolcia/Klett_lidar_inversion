%% ECSIM Lidar inversion - plots only
% close all
max_z = max(z_uv(id_ct_radar_uv));
nt = size(Ze_uv(:,:),1);
along_track = linspace(1,nt,nt);
%% Plot time series of measurements
TitleFigure=['ECSIM_time_series'];
figure('name', TitleFigure, 'NumberTitle','off', ...
    'units','centimeters','Position',[2 50 15 20]);
subplot(2,1,1)
Ze2= Ze_uv;
Ze2(Ze2 ==0)=NaN;
p11 = pcolor(along_track,z_uv,Ze2');
set(p11, 'EdgeColor', 'none');
set(gca,'Ydir','normal');
set(gca,'ylim',[0 max_z])
ylabel('z [m]')
hold on
l2 = plot(along_track,z_uv(norm_down_uv),'k','LineWidth',1.5);
hold on
l3 = plot(along_track,z_uv(id_cb_lidar_uv),'m','LineWidth',1.5);
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
p12 = pcolor(along_track,z_uv,Ba_MS_uv');
set(p12, 'EdgeColor', 'none');
set(gca,'Ydir','normal');
set(gca,'ylim',[0 max_z])
hold on
plot(along_track,z_uv(norm_down_uv),'k','LineWidth',1.5)
hold on
plot(along_track,z_uv(id_cb_lidar_uv),'m','LineWidth',1.5)
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
saveas(gcf,sprintf(TitleFigure),'epsc')
export_fig(sprintf(TitleFigure), '-eps', '-transparent')


%% Plot time series of extinction - the difference of the retrieved and original
TitleFigure=['Error_tau'];
figure('name', TitleFigure, 'NumberTitle','off', ...
    'units','centimeters','Position',[2 50 15 20]);
subplot(2,1,1)
tau_org3 = tau_uv;
tau_org3(tau_org3==0)=NaN;
p12 = pcolor(along_track,z_uv,tau_org3');
set(p12, 'EdgeColor', 'none');
set(gca,'Ydir','normal');
set(gca,'ylim',[0 600])
hold on
l1 = plot(along_track,z_uv(norm_down_uv),'k','LineWidth',1.5);
hold on
l2 = plot(along_track,z_uv(id_cb_lidar_uv),'m','LineWidth',1.5);
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

p11 = pcolor(along_track,z_uv,err_alpha_MS_corr_plot_uv');
set(p11, 'EdgeColor', 'none');
set(gca,'Ydir','normal');
set(gca,'ylim',[0 600])
hold on
plot(along_track,z_uv(norm_down_uv),'k','LineWidth',1.5)
hold on
plot(along_track,z_uv(id_cb_lidar_uv),'m','LineWidth',1.5)
ylabel('z [m]')
set(gca, 'FontSize',12)
c=colorbar('southoutside');
caxis([0 30])
c.Label.FontSize = 12;
c.Label.String = '%';
t = title('Retrieval Percent Error - with MS and RES correction', 'FontSize',12,'FontWeight','normal');
set(t, 'horizontalAlignment', 'left')
set(t, 'units', 'normalized')
legend([l1 l2],'Start of the normalisation interval', 'Cloud base estimate')
h1 = get(t, 'position');
set(t, 'position', [0 h1(2) h1(3)])
grid on
saveas(gcf,sprintf(TitleFigure),'epsc')
export_fig(sprintf(TitleFigure), '-eps', '-transparent') 

