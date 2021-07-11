% Lidar inversion - stats and tables 15m resolution
%% Calculate statistical data and make tables

for it = 1:length(id_cb_lidar_uv)
retrieved_cloud_uv{it} = id_cb_lidar_uv(it):norm_down_uv(it);
height_cloud_uv{it} = z(retrieved_cloud_uv{it});
end

for it = 1:size(alpha_ret_ms2ss_uv,1)
    for iih = 1:size(alpha_ret_ms2ss_uv,2)
      accu_alpha_plot_uv(it,iih) = (alpha_ret_ms2ss_uv(it,iih)./ Ext_uv(it,iih));
      err_alpha_MS_corr_plot_uv(it,iih) = (abs(alpha_ret_ms2ss_uv(it,iih) ...
     - Ext_uv(it,iih))./ Ext_uv(it,iih))*100; 
      err_alpha_MS_plot_uv(it,iih) = (abs(alpha_ret_ms2ss_nc_uv(it,iih) ...
     - Ext_uv(it,iih))./ Ext_uv(it,iih))*100; 
    
 accu_alpha_ss_uv(it,iih) = (alpha_ret_ss_uv(it,iih)./ Ext_uv(it,iih)) *100;
 err_alpha_ss_uv(it,iih) = (abs(alpha_ret_ss_uv(it,iih) ...
     - Ext_uv(it,iih))./ Ext_uv(it,iih))*100;
 
 accu_alpha_ss_nc_uv(it,iih) = (alpha_ret_ss_nc_uv(it,iih)./ Ext_uv(it,iih)) *100;
 err_alpha_ss_nc_uv(it,iih) = (abs(alpha_ret_ss_nc_uv(it,iih) ...
     - Ext_uv(it,iih))./ Ext_uv(it,iih))*100;
 
 accu_alpha_ms_nc_uv(it,iih) = (alpha_ret_ms_nc_uv(it,iih)./ Ext_uv(it,iih)) *100;
 err_alpha_ms_nc_uv(it,iih) = (abs(alpha_ret_ms_nc_uv(it,iih) ...
     - Ext_uv(it,iih))./ Ext_uv(it,iih))*100;
 
 accu_alpha_ms_uv(it,iih) = (alpha_ret_ms_uv(it,iih) ./ Ext_uv(it,iih)) *100;
 err_alpha_ms_uv(it,iih) = (abs(alpha_ret_ms_uv(it,iih) ...
     - Ext_uv(it,iih))./ Ext_uv(it,iih))*100;
 
 accu_alpha_ms_ss_uv(it,iih) = (alpha_ret_ms2ss_nc_uv(it,iih)./ Ext_uv(it,iih)) *100;
 err_alpha_ms_ss_uv(it,iih) = (abs(alpha_ret_ms2ss_nc_uv(it,iih) ...
     - Ext_uv(it,iih))./ Ext_uv(it,iih))*100;
   
 accu_alpha_ms_ss_res_uv(it,iih) = (alpha_ret_ms2ss_uv(it,iih)./ Ext_uv(it,iih)) *100;
 err_alpha_ms_ss_res_uv(it,iih) = (abs(alpha_ret_ms2ss_uv(it,iih) ...
    - Ext_uv(it,iih))./ Ext_uv(it,iih))*100;

     end  

accu_alpha_plot_uv(it,1:retrieved_cloud_uv{it}-1)=NaN;
accu_alpha_plot_uv(it,max(retrieved_cloud_uv{it}):end)=NaN;
err_alpha_MS_corr_plot_uv(it,1:retrieved_cloud_uv{it}-1)=NaN;
err_alpha_MS_corr_plot_uv(it,max(retrieved_cloud_uv{it}):end)=NaN;
err_alpha_MS_plot_uv(it,1:retrieved_cloud_uv{it}-1)=NaN;
err_alpha_MS_plot_uv(it,max(retrieved_cloud_uv{it}):end)=NaN;

accu_alpha_ss_uv(it,1:retrieved_cloud_uv{it}-1)=NaN;
accu_alpha_ss_uv(it,max(retrieved_cloud_uv{it}):end)=NaN;
err_alpha_ss_uv(it,1:retrieved_cloud_uv{it}-1)=NaN;
err_alpha_ss_uv(it,max(retrieved_cloud_uv{it}):end)=NaN;

accu_alpha_ss_nc_uv(it,1:retrieved_cloud_uv{it}-1)=NaN;
accu_alpha_ss_nc_uv(it,max(retrieved_cloud_uv{it}):end)=NaN;
err_alpha_ss_nc_uv(it,1:retrieved_cloud_uv{it}-1)=NaN;
err_alpha_ss_nc_uv(it,max(retrieved_cloud_uv{it}):end)=NaN;

accu_alpha_ms_uv(it,1:retrieved_cloud_uv{it}-1)=NaN;
accu_alpha_ms_uv(it,max(retrieved_cloud_uv{it}):end)=NaN;
err_alpha_ms_uv(it,1:retrieved_cloud_uv{it}-1)=NaN;
err_alpha_ms_uv(it,max(retrieved_cloud_uv{it}):end)=NaN;

accu_alpha_ms_ss_uv(it,1:retrieved_cloud_uv{it}-1)=NaN;
accu_alpha_ms_ss_uv(it,max(retrieved_cloud_uv{it}):end)=NaN;
err_alpha_ms_ss_uv(it,1:retrieved_cloud_uv{it}-1)=NaN;
err_alpha_ms_ss_uv(it,max(retrieved_cloud_uv{it}):end)=NaN;

accu_alpha_ms_ss_res_uv(it,1:retrieved_cloud_uv{it}-1)=NaN;
accu_alpha_ms_ss_res_uv(it,max(retrieved_cloud_uv{it}):end)=NaN;
err_alpha_ms_ss_res_uv(it,1:retrieved_cloud_uv{it}-1)=NaN;
err_alpha_ms_ss_res_uv(it,max(retrieved_cloud_uv{it}):end)=NaN;

end

accu_alpha_ss_uv(accu_alpha_ss_uv ==0) = NaN;
err_alpha_ss_uv(err_alpha_ss_uv==0)= NaN;

accu_alpha_ss_nc_uv(accu_alpha_ss_nc_uv ==0) = NaN;
err_alpha_ss_nc_uv(err_alpha_ss_nc_uv==0)= NaN;

accu_alpha_ms_nc_uv( accu_alpha_ms_nc_uv ==0) = NaN;
err_alpha_ms_nc_uv(err_alpha_ms_nc_uv==0)= NaN;

accu_alpha_ms_uv( accu_alpha_ms_uv ==0) = NaN;
err_alpha_ms_uv(err_alpha_ms_uv==0)= NaN;

accu_alpha_ms_ss_uv(accu_alpha_ms_ss_uv==0) = NaN;
err_alpha_ms_ss_uv(err_alpha_ms_ss_uv==0)= NaN;

accu_alpha_ms_ss_res_uv(accu_alpha_ms_ss_res_uv==0) = NaN;
err_alpha_ms_ss_res_uv(err_alpha_ms_ss_res_uv==0)= NaN;

accu_alpha_ss_uv(accu_alpha_ss_uv == Inf) = NaN;
err_alpha_ss_uv(err_alpha_ss_uv == Inf)= NaN;

accu_alpha_ss_nc_uv(accu_alpha_ss_nc_uv ==Inf) = NaN;
err_alpha_ss_nc_uv(err_alpha_ss_nc_uv==Inf)= NaN;

accu_alpha_ms_nc_uv( accu_alpha_ms_nc_uv ==Inf) = NaN;
err_alpha_ms_nc_uv(err_alpha_ms_nc_uv==Inf)= NaN;

accu_alpha_ms_uv( accu_alpha_ms_uv ==Inf) = NaN;
err_alpha_ms_uv(err_alpha_ms_uv==Inf)= NaN;

accu_alpha_ms_ss_uv(accu_alpha_ms_ss_uv==Inf) = NaN;
err_alpha_ms_ss_uv(err_alpha_ms_ss_uv==Inf)= NaN;

accu_alpha_ms_ss_res_uv(accu_alpha_ms_ss_res_uv==Inf) = NaN;
err_alpha_ms_ss_res_uv(err_alpha_ms_ss_res_uv==Inf)= NaN;


%% Latex table
% fprintf('m from cb & ACU   & ERR  -- ECSIM SS retrieved \\ \n')
% for ii=1:length(retrieved_cloud_uv{1})
%     fprintf('%8.1f & %3.2f%% &  %3.2f%%  \\\\ \n', ...
%         z(retrieved_cloud_uv{1}(ii)) - z(retrieved_cloud_uv{1}(1)) , ...
%         nanmean(accu_alpha_ss_nc(:,ii)), ...
%         nanmean(err_alpha_ss_nc(:,ii)))
% end
% % 
% % fprintf('m from cb & ACU   & ERR  -- ECSIM SS retrieved with RES correction \\ \n')
% % for ii=1:length(retrieved_cloud_uv{1})
% %     fprintf('%8.1f & %3.2f%% &  %3.2f%%  \\\\ \n', ...
% %         z(retrieved_cloud_uv{1}(ii)) - z(retrieved_cloud_uv{1}(1)) , ...
% %         nanmean(accu_alpha_ss(:,ii)), ...
% %         nanmean(err_alpha_ss(:,ii)))
% % end
% % 
% % fprintf('m from cb & ACU   & ERR  -- ECSIM MS \\ \n')
% % for ii=1:length(retrieved_cloud_uv{1})
% %     fprintf('%8.1f & %3.2f%% &  %3.2f%%  \\\\ \n', ...
% %         z(retrieved_cloud_uv{1}(ii)) - z(retrieved_cloud_uv{1}(1)) , ...
% %         nanmean(accu_alpha_ms_nc(:,ii)), ...
% %         nanmean(err_alpha_ms_nc(:,ii)))
% % end
% % 
% % fprintf('m from cb & ACU   & ERR  -- ECSIM MS with RES correction \\ \n')
% % for ii=1:length(retrieved_cloud_uv{1})
% %     fprintf('%8.1f & %3.2f%% &  %3.2f%%  \\\\ \n', ...
% %         z(retrieved_cloud_uv{1}(ii)) - z(retrieved_cloud_uv{1}(1)) , ...
% %         nanmean(accu_alpha_ms(:,ii)), ...
% %         nanmean(err_alpha_ms(:,ii)))
% % end
% % 
% % fprintf('m from cb & ACU   & ERR  -- ECSIM MS with SS correction \\ \n')
% % for ii=1:length(retrieved_cloud_uv{1})
% %     fprintf('%8.1f & %3.2f%% &  %3.2f%%  \\\\ \n', ...
% %         z(retrieved_cloud_uv{1}(ii)) - z(retrieved_cloud_uv{1}(1)) , ...
% %         nanmean(accu_alpha_ms_ss(:,ii)), ...
% %         nanmean(err_alpha_ms_ss(:,ii)))
% % end
% % 
% % fprintf('m from cb & ACU   & ERR  -- ECSIM MS with SS & RES correction \\ \n')
% % for ii=1:length(retrieved_cloud_uv{1})
% %     fprintf('%8.1f & %3.2f%% &  %3.2f%%  \\\\ \n', ...
% %         z(retrieved_cloud_uv{1}(ii)) - z(retrieved_cloud_uv{1}(1)) , ...
% %         nanmean(accu_alpha_ms_ss_res(:,ii)), ...
% %         nanmean(err_alpha_ms_ss_res(:,ii)))
% % end
% fprintf('m from cb & ACU   & ERR  -- ECSIM SS retrieved \\ \n')
% for ii=1:length(retrieved_cloud_uv{1})
%     fprintf('%8.1f & %3.2f%% &  %3.2f%%  \\\\ \n', ...
%         z(retrieved_cloud_uv{1}(ii)) - z(retrieved_cloud_uv{1}(1)) , ...
%         nanmean(accu_alpha_ss_nc(:,retrieved_cloud_uv{1}(ii))), ...
%         nanmean(err_alpha_ss_nc(:,retrieved_cloud_uv{1}(ii))))
% end
% 
% for ii=1:length(retrieved_cloud_uv{1})
%     fprintf('%8.1f & %3.2f%% &  %3.2f%%  \\\\ \n', ...
%         z(retrieved_cloud_uv{1}(ii)) - z(retrieved_cloud_uv{1}(1)) , ...
%         nanmean(accu_alpha_ss_nc(:,retrieved_cloud_uv{1}(ii))), ...
%         nanmean(err_alpha_ss_nc(:,retrieved_cloud_uv{1}(ii))))
% end
% 
% fprintf('m from cb & ACU   & ERR  -- ECSIM SS retrieved with RES correction \\ \n')
% for ii=1:length(retrieved_cloud_uv{1})
%     fprintf('%8.1f & %3.2f%% &  %3.2f%%  \\\\ \n', ...
%         z(retrieved_cloud_uv{1}(ii)) - z(retrieved_cloud_uv{1}(1)) , ...
%         nanmean(accu_alpha_ss(:,ii)), ...
%         nanmean(err_alpha_ss(:,ii)))
% end
% 
% fprintf('m from cb & ACU   & ERR  -- ECSIM MS \\ \n')
% for ii=1:length(retrieved_cloud_uv{1})
%     fprintf('%8.1f & %3.2f%% &  %3.2f%%  \\\\ \n', ...
%         z(retrieved_cloud_uv{1}(ii)) - z(retrieved_cloud_uv{1}(1)) , ...
%         nanmean(accu_alpha_ms_nc(:,retrieved_cloud_uv{1}(ii))), ...
%         nanmean(err_alpha_ms_nc(:,retrieved_cloud_uv{1}(ii))))
% end
% 
% fprintf('m from cb & ACU   & ERR  -- ECSIM MS with RES correction \\ \n')
% for ii=1:length(retrieved_cloud_uv{1})
%     fprintf('%8.1f & %3.2f%% &  %3.2f%%  \\\\ \n', ...
%         z(retrieved_cloud_uv{1}(ii)) - z(retrieved_cloud_uv{1}(1)) , ...
%         nanmean(accu_alpha_ms(:,retrieved_cloud_uv{1}(ii))), ...
%         nanmean(err_alpha_ms(:,retrieved_cloud_uv{1}(ii))))
% end
% 
% fprintf('m from cb & ACU   & ERR  -- ECSIM MS with SS correction \\ \n')
% for ii=1:length(retrieved_cloud_uv{1})
%     fprintf('%8.1f & %3.2f%% &  %3.2f%%  \\\\ \n', ...
%         z(retrieved_cloud_uv{1}(ii)) - z(retrieved_cloud_uv{1}(1)) , ...
%         nanmean(accu_alpha_ms_ss(:,retrieved_cloud_uv{1}(ii))), ...
%         nanmean(err_alpha_ms_ss(:,retrieved_cloud_uv{1}(ii))))
% end
% 
fprintf('m from cb & ACU   & ERR  -- ECSIM MS with SS correction \\ \n')
for ii=1:length(retrieved_cloud_uv{1})
    fprintf('%8.1f & %3.2f%% &  %3.2f%%  \\\\ \n', ...
        z_uv(retrieved_cloud_uv{1}(ii)) - z_uv(retrieved_cloud_uv{1}(1)) , ...
        nanmean(accu_alpha_ms_ss_uv(:,retrieved_cloud_uv{1}(ii))), ...
        nanmean(err_alpha_ms_ss_uv(:,retrieved_cloud_uv{1}(ii))))
end

fprintf('m from cb & ACU   & ERR  -- ECSIM MS with SS & RES correction \\ \n')
for ii=1:length(retrieved_cloud_uv{1})
    fprintf('%8.1f & %3.2f%% &  %3.2f%%  \\\\ \n', ...
        z_uv(retrieved_cloud_uv{1}(ii)) - z_uv(retrieved_cloud_uv{1}(1)) , ...
        nanmean(accu_alpha_ms_ss_res_uv(:,retrieved_cloud_uv{1}(ii))), ...
        nanmean(err_alpha_ms_ss_res_uv(:,retrieved_cloud_uv{1}(ii))))
end

