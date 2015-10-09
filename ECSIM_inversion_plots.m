%% ECSIM Lidar inversion - plots only
max_height = max(height(norm_up + 15));

%% Plot time series of measurements
TitleFigure=['ECSIM data time series'];
figure('name', TitleFigure, 'NumberTitle','off', ...
    'units','centimeters','Position',[2 50 15 20]);
subplot(2,1,1)
p11 = pcolor(time,height,Ze');
set(p11, 'EdgeColor', 'none');
set(gca,'Ydir','normal');
set(gca,'ylim',[0 max_height])
ylabel('Height [m]')
hold on
l1 = plot(time,height(norm_down),'r');
% set(gca,'YTick',500:1000:max_height)
% set(gca, 'YTickLabel', [0.5 1.5 2.5]); 
set(gca, 'FontSize',8)
c=colorbar('southoutside');
c.Label.String = 'Z [mm^6/m^3]';
c.Label.FontSize = 8;
% xlabel('Time [UTC]')
t = title('Radar Reflectivity Factor', 'FontSize',10,'FontWeight','normal');
set(t, 'horizontalAlignment', 'left')
set(t, 'units', 'normalized')
h1 = get(t, 'position');
set(t, 'position', [0 h1(2) h1(3)])
grid on
legend(l1, 'begining of the normalisation interval')

% ax = gca;
% axpos = ax.Position;
% cpos = c.Position;
% cpos(1) = axpos(1);
% cpos(2) = cpos(2) + 0.045;
% cpos(3) = axpos(3);
% cpos(4) = 0.5*cpos(4);
% c.Position = cpos;
% ax.Position = axpos;
% clear ax cpos c t

subplot(2,1,2)
p12 = pcolor(time,height,beta_atten');
set(p12, 'EdgeColor', 'none');
set(gca,'Ydir','normal');
set(gca,'ylim',[0 max_height])
hold on
l2 = plot(time,height(norm_down),'r');
hold on
l3= plot(time,height(id_cb_lidar),'g');
hold on
l4 = plot(time,height(id_cb_radar),'m');
% xlabel('Time [UTC]')
ylabel('Height [m]')
% set(gca,'YTick',500:1000:max_height)
% set(gca, 'YTickLabel', [1 1.5 2.5]); 
set(gca, 'FontSize',8)
t = title('Attenuated Backscatter Coefficient', 'FontSize',10,'FontWeight','normal');
set(t, 'horizontalAlignment', 'left')
set(t, 'units', 'normalized')
h1 = get(t, 'position');
set(t, 'position', [0 h1(2) h1(3)])
legend([l2,l3,l4], 'begining of the normalisation interval', 'cloudbase from lidar', ...
    'cloud base from radar')
% caxis([0 5*1.e-5])
c=colorbar('southoutside');
c.Label.String = '\beta [m^{-1} sr^{-1}]';
% % set(c,'YTick',1*1.e-5:1*1.e-5:4*1.e-5)
% set(c, 'YTickLabel', cellstr(num2str(reshape(get(c, 'YTick'),[],1),'%2.e')) )
% c.Label.FontSize = 8;
% grid on
% ax = gca;
% axpos = ax.Position;
% cpos = c.Position;
% cpos(1) = axpos(1);
% cpos(2) = cpos(2) + 0.045;
% cpos(3) = axpos(3);
% cpos(4) = 0.5*cpos(4);
% c.Position = cpos;
% ax.Position = axpos;
% clear ax cpos c t


%% Plot time series of extinction
TitleFigure=['Original and Retrieved Extinction'];
figure('name', TitleFigure, 'NumberTitle','off', ...
    'units','centimeters','Position',[2 50 15 20]);
subplot(2,1,1)
p11 = pcolor(time,height,Ext');
set(p11, 'EdgeColor', 'none');
set(gca,'Ydir','normal');
set(gca,'ylim',[0 max_height])
hold on
plot(time,height(norm_down),'r')
hold on
plot(time,height(id_cb_lidar),'g')
ylabel('Height [m]')
% set(gca,'YTick',500:1000:max_height)
% set(gca, 'YTickLabel', [0.5 1.5 2.5]); 
set(gca, 'FontSize',8)
c=colorbar('southoutside');
caxis([0 0.06])
c.Label.String = '[m^{-1}]';
c.Label.FontSize = 8;
% xlabel('Time [UTC]')
t = title('Modeled Extinction from ECSIM', 'FontSize',10,'FontWeight','normal');
set(t, 'horizontalAlignment', 'left')
set(t, 'units', 'normalized')
h1 = get(t, 'position');
set(t, 'position', [0 h1(2) h1(3)])
grid on
% ax = gca;
% axpos = ax.Position;
% cpos = c.Position;
% cpos(1) = axpos(1);
% cpos(2) = cpos(2) + 0.045;
% cpos(3) = axpos(3);
% cpos(4) = 0.5*cpos(4);
% c.Position = cpos;
% ax.Position = axpos;
% clear ax cpos c t

subplot(2,1,2)
p12 = pcolor(time,height,alpha');
set(p12, 'EdgeColor', 'none');
set(gca,'Ydir','normal');
set(gca,'ylim',[0 max_height])
% hold on
% plot(time,height(norm_down),'r')
% hold on
% plot(time,height(id_cb_lidar),'g')
% xlabel('Time [UTC]')
ylabel('Height [m]')
% set(gca,'YTick',500:1000:max_height)
% set(gca, 'YTickLabel', [1 1.5 2.5]); 
set(gca, 'FontSize',8)
t = title('Extinction Retrieved', 'FontSize',10,'FontWeight','normal');
set(t, 'horizontalAlignment', 'left')
set(t, 'units', 'normalized')
h1 = get(t, 'position');
set(t, 'position', [0 h1(2) h1(3)])
caxis([0 0.06])
c=colorbar('southoutside');
c.Label.String = '[m^{-1}]';
% 
% Plot time series of extinction - the difference of the retrieved and original
TitleFigure=['Original and Retrieved Extinction'];
figure('name', TitleFigure, 'NumberTitle','off', ...
    'units','centimeters','Position',[2 50 15 10]);
subplot(1,1,1)
p11 = pcolor(time,height,accu_alpha_plot');
set(p11, 'EdgeColor', 'none');
set(gca,'Ydir','normal');
set(gca,'ylim',[0 max_height])
hold on
plot(time,height(norm_down),'r')
hold on
plot(time,height(id_cb_lidar),'g')
ylabel('Height [m]')
set(gca,'YTick',500:1000:max_height)
set(gca, 'YTickLabel', [0.5 1.5 2.5]); 
set(gca, 'FontSize',8)
c=colorbar('southoutside');
caxis([0 1.1])
c.Label.FontSize = 8;
xlabel('Time [UTC]')
t = title('Retrieved extinction divided by original', 'FontSize',10,'FontWeight','normal');
set(t, 'horizontalAlignment', 'left')
set(t, 'units', 'normalized')
h1 = get(t, 'position');
set(t, 'position', [0 h1(2) h1(3)])
grid on
% ax = gca;
% axpos = ax.Position;
% cpos = c.Position;
% cpos(1) = axpos(1);
% cpos(2) = cpos(2) + 0.65;
% cpos(3) = axpos(3);
% cpos(4) = 0.5*cpos(4);
% c.Position = cpos;
% ax.Position = axpos;
% clear ax cpos c t


% subplot(2,1,2)
% p12 = pcolor(time,height_cloud{1},accu_alpha_sscat');
% p11 = pcolor(time,height(id_cb_lidar:norm_down), ...
%     (alpha(:,id_cb_lidar:norm_down)./Ext(:,id_cb_lidar:norm_down))');
% set(p11, 'EdgeColor', 'none');
% set(gca,'Ydir','normal');
% set(gca,'ylim',[0 max_height])
% hold on
% plot(time,height(norm_down),'r')
% hold on
% plot(time,height(id_cb_lidar),'g')
% ylabel('Height [m]')
% set(gca,'YTick',500:1000:max_height)
% set(gca, 'YTickLabel', [0.5 1.5 2.5]); 
% set(gca, 'FontSize',8)
% c=colorbar('southoutside');
% caxis([0 1])
% c.Label.String = '[m^{-1}]';
% c.Label.FontSize = 8;
% xlabel('Time [UTC]')
% t = title('Single scattering extinction divided by original extinction', 'FontSize',10,'FontWeight','normal');
% set(t, 'horizontalAlignment', 'left')
% set(t, 'units', 'normalized')
% h1 = get(t, 'position');
% set(t, 'position', [0 h1(2) h1(3)])
% grid on
% ax = gca;
% axpos = ax.Position;
% cpos = c.Position;
% cpos(1) = axpos(1);
% cpos(2) = cpos(2) + 0.045;
% cpos(3) = axpos(3);
% cpos(4) = 0.5*cpos(4);
% c.Position = cpos;
% ax.Position = axpos;
% clear ax cpos c t
% 
% % subplot(3,1,3)
% % p12 = pcolor(time,height(id_cb_lidar:norm_up),(alpha(:,id_cb_lidar:norm_up) - alpha_s(:,id_cb_lidar:norm_up))');
% % set(p12, 'EdgeColor', 'none');
% % set(gca,'Ydir','normal');
% % set(gca,'ylim',[0 max_height])
% % hold on
% % plot(time,height(norm_down),'r')
% % % xlabel('Time [UTC]')
% % ylabel('Height [m]')
% % % set(gca,'YTick',500:1000:max_height)
% % % set(gca, 'YTickLabel', [1 1.5 2.5]);
% % caxis([-0.0025 0.0025])
% % set(gca, 'FontSize',8)
% % t = title('Retrieve Extinction witout MS minus Retrieved Extinction with MS', 'FontSize',10,'FontWeight','normal');
% % set(t, 'horizontalAlignment', 'left')
% % set(t, 'units', 'normalized')
% % h1 = get(t, 'position');
% % set(t, 'position', [0 h1(2) h1(3)])
% % % caxis([-5 5])
% % c=colorbar('southoutside');
% % c.Label.String = '[m^{-1}]';


