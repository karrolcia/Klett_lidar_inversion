% Lidar inversion - stats and tables 2.5m resolution
%% Calculate statistical data and make tables

for it = 1:length(id_cb_lidar)
retrieved_cloud{it} = id_cb_lidar(it):norm_down(it);
height_cloud{it} = z(retrieved_cloud{it});
end

for it = 1:size(alpha_ret_ms2ss,1)
    for iih = 1:size(alpha_ret_ms2ss,2)
      accu_alpha_plot(it,iih) = (alpha_ret_ms2ss(it,iih)./ Ext(it,iih));
      err_alpha_MS_corr_plot(it,iih) = (abs(alpha_ret_ms2ss(it,iih) ...
     - Ext(it,iih))./ Ext(it,iih))*100; 
      err_alpha_MS_plot(it,iih) = (abs(alpha_ret_ms2ss_nc(it,iih) ...
     - Ext(it,iih))./ Ext(it,iih))*100; 
    
 accu_alpha_ss(it,iih) = (alpha_ret_ss(it,iih)./ Ext(it,iih)) *100;
 err_alpha_ss(it,iih) = (abs(alpha_ret_ss(it,iih) ...
     - Ext(it,iih))./ Ext(it,iih))*100;
 
 accu_alpha_ss_nc(it,iih) = (alpha_ret_ss_nc(it,iih)./ Ext(it,iih)) *100;
 err_alpha_ss_nc(it,iih) = (abs(alpha_ret_ss_nc(it,iih) ...
     - Ext(it,iih))./ Ext(it,iih))*100;
 
 accu_alpha_ms_nc(it,iih) = (alpha_ret_ms_nc(it,iih)./ Ext(it,iih)) *100;
 err_alpha_ms_nc(it,iih) = (abs(alpha_ret_ms_nc(it,iih) ...
     - Ext(it,iih))./ Ext(it,iih))*100;
 
 accu_alpha_ms(it,iih) = (alpha_ret_ms(it,iih) ./ Ext(it,iih)) *100;
 err_alpha_ms(it,iih) = (abs(alpha_ret_ms(it,iih) ...
     - Ext(it,iih))./ Ext(it,iih))*100;
 
 accu_alpha_ms_ss(it,iih) = (alpha_ret_ms2ss_nc(it,iih)./ Ext(it,iih)) *100;
 err_alpha_ms_ss(it,iih) = (abs(alpha_ret_ms2ss_nc(it,iih) ...
     - Ext(it,iih))./ Ext(it,iih))*100;
   
 accu_alpha_ms_ss_res(it,iih) = (alpha_ret_ms2ss(it,iih)./ Ext(it,iih)) *100;
 err_alpha_ms_ss_res(it,iih) = (abs(alpha_ret_ms2ss(it,iih) ...
    - Ext(it,iih))./ Ext(it,iih))*100;

     end  

accu_alpha_plot(it,1:retrieved_cloud{it}-1)=NaN;
accu_alpha_plot(it,max(retrieved_cloud{it}):end)=NaN;
err_alpha_MS_corr_plot(it,1:retrieved_cloud{it}-1)=NaN;
err_alpha_MS_corr_plot(it,max(retrieved_cloud{it}):end)=NaN;
err_alpha_MS_plot(it,1:retrieved_cloud{it}-1)=NaN;
err_alpha_MS_plot(it,max(retrieved_cloud{it}):end)=NaN;

accu_alpha_ss(it,1:retrieved_cloud{it}-1)=NaN;
accu_alpha_ss(it,max(retrieved_cloud{it}):end)=NaN;
err_alpha_ss(it,1:retrieved_cloud{it}-1)=NaN;
err_alpha_ss(it,max(retrieved_cloud{it}):end)=NaN;

accu_alpha_ss_nc(it,1:retrieved_cloud{it}-1)=NaN;
accu_alpha_ss_nc(it,max(retrieved_cloud{it}):end)=NaN;
err_alpha_ss_nc(it,1:retrieved_cloud{it}-1)=NaN;
err_alpha_ss_nc(it,max(retrieved_cloud{it}):end)=NaN;

accu_alpha_ms(it,1:retrieved_cloud{it}-1)=NaN;
accu_alpha_ms(it,max(retrieved_cloud{it}):end)=NaN;
err_alpha_ms(it,1:retrieved_cloud{it}-1)=NaN;
err_alpha_ms(it,max(retrieved_cloud{it}):end)=NaN;

accu_alpha_ms_ss(it,1:retrieved_cloud{it}-1)=NaN;
accu_alpha_ms_ss(it,max(retrieved_cloud{it}):end)=NaN;
err_alpha_ms_ss(it,1:retrieved_cloud{it}-1)=NaN;
err_alpha_ms_ss(it,max(retrieved_cloud{it}):end)=NaN;

accu_alpha_ms_ss_res(it,1:retrieved_cloud{it}-1)=NaN;
accu_alpha_ms_ss_res(it,max(retrieved_cloud{it}):end)=NaN;
err_alpha_ms_ss_res(it,1:retrieved_cloud{it}-1)=NaN;
err_alpha_ms_ss_res(it,max(retrieved_cloud{it}):end)=NaN;

end

accu_alpha_ss(accu_alpha_ss ==0) = NaN;
err_alpha_ss(err_alpha_ss==0)= NaN;

accu_alpha_ss_nc(accu_alpha_ss_nc ==0) = NaN;
err_alpha_ss_nc(err_alpha_ss_nc==0)= NaN;

accu_alpha_ms_nc( accu_alpha_ms_nc ==0) = NaN;
err_alpha_ms_nc(err_alpha_ms_nc==0)= NaN;

accu_alpha_ms( accu_alpha_ms ==0) = NaN;
err_alpha_ms(err_alpha_ms==0)= NaN;

accu_alpha_ms_ss(accu_alpha_ms_ss==0) = NaN;
err_alpha_ms_ss(err_alpha_ms_ss==0)= NaN;

accu_alpha_ms_ss_res(accu_alpha_ms_ss_res==0) = NaN;
err_alpha_ms_ss_res(err_alpha_ms_ss_res==0)= NaN;

accu_alpha_ss(accu_alpha_ss == Inf) = NaN;
err_alpha_ss(err_alpha_ss == Inf)= NaN;

accu_alpha_ss_nc(accu_alpha_ss_nc ==Inf) = NaN;
err_alpha_ss_nc(err_alpha_ss_nc==Inf)= NaN;

accu_alpha_ms_nc( accu_alpha_ms_nc ==Inf) = NaN;
err_alpha_ms_nc(err_alpha_ms_nc==Inf)= NaN;

accu_alpha_ms( accu_alpha_ms ==Inf) = NaN;
err_alpha_ms(err_alpha_ms==Inf)= NaN;

accu_alpha_ms_ss(accu_alpha_ms_ss==Inf) = NaN;
err_alpha_ms_ss(err_alpha_ms_ss==Inf)= NaN;

accu_alpha_ms_ss_res(accu_alpha_ms_ss_res==Inf) = NaN;
err_alpha_ms_ss_res(err_alpha_ms_ss_res==Inf)= NaN;


%% Latex table
% fprintf('m from cb & ACU   & ERR  -- ECSIM SS retrieved \\ \n')
% for ii=1:length(retrieved_cloud{1})
%     fprintf('%8.1f & %3.2f%% &  %3.2f%%  \\\\ \n', ...
%         z(retrieved_cloud{1}(ii)) - z(retrieved_cloud{1}(1)) , ...
%         nanmean(accu_alpha_ss_nc(:,ii)), ...
%         nanmean(err_alpha_ss_nc(:,ii)))
% end
% % 
% % fprintf('m from cb & ACU   & ERR  -- ECSIM SS retrieved with RES correction \\ \n')
% % for ii=1:length(retrieved_cloud{1})
% %     fprintf('%8.1f & %3.2f%% &  %3.2f%%  \\\\ \n', ...
% %         z(retrieved_cloud{1}(ii)) - z(retrieved_cloud{1}(1)) , ...
% %         nanmean(accu_alpha_ss(:,ii)), ...
% %         nanmean(err_alpha_ss(:,ii)))
% % end
% % 
% % fprintf('m from cb & ACU   & ERR  -- ECSIM MS \\ \n')
% % for ii=1:length(retrieved_cloud{1})
% %     fprintf('%8.1f & %3.2f%% &  %3.2f%%  \\\\ \n', ...
% %         z(retrieved_cloud{1}(ii)) - z(retrieved_cloud{1}(1)) , ...
% %         nanmean(accu_alpha_ms_nc(:,ii)), ...
% %         nanmean(err_alpha_ms_nc(:,ii)))
% % end
% % 
% % fprintf('m from cb & ACU   & ERR  -- ECSIM MS with RES correction \\ \n')
% % for ii=1:length(retrieved_cloud{1})
% %     fprintf('%8.1f & %3.2f%% &  %3.2f%%  \\\\ \n', ...
% %         z(retrieved_cloud{1}(ii)) - z(retrieved_cloud{1}(1)) , ...
% %         nanmean(accu_alpha_ms(:,ii)), ...
% %         nanmean(err_alpha_ms(:,ii)))
% % end
% % 
% % fprintf('m from cb & ACU   & ERR  -- ECSIM MS with SS correction \\ \n')
% % for ii=1:length(retrieved_cloud{1})
% %     fprintf('%8.1f & %3.2f%% &  %3.2f%%  \\\\ \n', ...
% %         z(retrieved_cloud{1}(ii)) - z(retrieved_cloud{1}(1)) , ...
% %         nanmean(accu_alpha_ms_ss(:,ii)), ...
% %         nanmean(err_alpha_ms_ss(:,ii)))
% % end
% % 
% % fprintf('m from cb & ACU   & ERR  -- ECSIM MS with SS & RES correction \\ \n')
% % for ii=1:length(retrieved_cloud{1})
% %     fprintf('%8.1f & %3.2f%% &  %3.2f%%  \\\\ \n', ...
% %         z(retrieved_cloud{1}(ii)) - z(retrieved_cloud{1}(1)) , ...
% %         nanmean(accu_alpha_ms_ss_res(:,ii)), ...
% %         nanmean(err_alpha_ms_ss_res(:,ii)))
% % end
% fprintf('m from cb & ACU   & ERR  -- ECSIM SS retrieved \\ \n')
% for ii=1:length(retrieved_cloud{1})
%     fprintf('%8.1f & %3.2f%% &  %3.2f%%  \\\\ \n', ...
%         z(retrieved_cloud{1}(ii)) - z(retrieved_cloud{1}(1)) , ...
%         nanmean(accu_alpha_ss_nc(:,retrieved_cloud{1}(ii))), ...
%         nanmean(err_alpha_ss_nc(:,retrieved_cloud{1}(ii))))
% end
% 
% for ii=1:length(retrieved_cloud{1})
%     fprintf('%8.1f & %3.2f%% &  %3.2f%%  \\\\ \n', ...
%         z(retrieved_cloud{1}(ii)) - z(retrieved_cloud{1}(1)) , ...
%         nanmean(accu_alpha_ss_nc(:,retrieved_cloud{1}(ii))), ...
%         nanmean(err_alpha_ss_nc(:,retrieved_cloud{1}(ii))))
% end
% 
% fprintf('m from cb & ACU   & ERR  -- ECSIM SS retrieved with RES correction \\ \n')
% for ii=1:length(retrieved_cloud{1})
%     fprintf('%8.1f & %3.2f%% &  %3.2f%%  \\\\ \n', ...
%         z(retrieved_cloud{1}(ii)) - z(retrieved_cloud{1}(1)) , ...
%         nanmean(accu_alpha_ss(:,ii)), ...
%         nanmean(err_alpha_ss(:,ii)))
% end
% 
% fprintf('m from cb & ACU   & ERR  -- ECSIM MS \\ \n')
% for ii=1:length(retrieved_cloud{1})
%     fprintf('%8.1f & %3.2f%% &  %3.2f%%  \\\\ \n', ...
%         z(retrieved_cloud{1}(ii)) - z(retrieved_cloud{1}(1)) , ...
%         nanmean(accu_alpha_ms_nc(:,retrieved_cloud{1}(ii))), ...
%         nanmean(err_alpha_ms_nc(:,retrieved_cloud{1}(ii))))
% end
% 
% fprintf('m from cb & ACU   & ERR  -- ECSIM MS with RES correction \\ \n')
% for ii=1:length(retrieved_cloud{1})
%     fprintf('%8.1f & %3.2f%% &  %3.2f%%  \\\\ \n', ...
%         z(retrieved_cloud{1}(ii)) - z(retrieved_cloud{1}(1)) , ...
%         nanmean(accu_alpha_ms(:,retrieved_cloud{1}(ii))), ...
%         nanmean(err_alpha_ms(:,retrieved_cloud{1}(ii))))
% end
% 
% fprintf('m from cb & ACU   & ERR  -- ECSIM MS with SS correction \\ \n')
% for ii=1:length(retrieved_cloud{1})
%     fprintf('%8.1f & %3.2f%% &  %3.2f%%  \\\\ \n', ...
%         z(retrieved_cloud{1}(ii)) - z(retrieved_cloud{1}(1)) , ...
%         nanmean(accu_alpha_ms_ss(:,retrieved_cloud{1}(ii))), ...
%         nanmean(err_alpha_ms_ss(:,retrieved_cloud{1}(ii))))
% end
% 
fprintf('m from cb & ACU   & ERR  -- ECSIM MS with SS correction \\ \n')
for ii=1:length(retrieved_cloud{1})
    fprintf('%8.1f & %3.2f%% &  %3.2f%%  \\\\ \n', ...
        z(retrieved_cloud{1}(ii)) - z(retrieved_cloud{1}(1)) , ...
        nanmean(accu_alpha_ms_ss(:,retrieved_cloud{1}(ii))), ...
        nanmean(err_alpha_ms_ss(:,retrieved_cloud{1}(ii))))
end

fprintf('m from cb & ACU   & ERR  -- ECSIM MS with SS & RES correction \\ \n')
for ii=1:length(retrieved_cloud{1})
    fprintf('%8.1f & %3.2f%% &  %3.2f%%  \\\\ \n', ...
        z(retrieved_cloud{1}(ii)) - z(retrieved_cloud{1}(1)) , ...
        nanmean(accu_alpha_ms_ss_res(:,retrieved_cloud{1}(ii))), ...
        nanmean(err_alpha_ms_ss_res(:,retrieved_cloud{1}(ii))))
end

