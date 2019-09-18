function []=fcn_plot_simul_anneal(data_init_optim,T_loss,nodes,sel_nodes,vars_show,thres_ind,plot_settings)


label_ticks_fontsize=plot_settings(1); label_fontsize=plot_settings(2);


if size(T_loss,2)>1 
    error_data=T_loss(1:thres_ind,vars_show);
    init_error=T_loss(1,2); best_error=T_loss(end,2);
else
    error_data=T_loss; 
    if isempty(vars_show); vars_show=2;end
end

if isempty(thres_ind)
    thres_ind=max(size(T_loss)); init_error=T_loss(1); best_error=T_loss(end);
end

fig_subpl1=subplot(1,2,1); 
plot(1:thres_ind, error_data,'LineWidth',4); xlim([1 thres_ind]); if init_error/best_error>30; set(gca,'yscale','log'); end
set(gca,'FontSize',label_ticks_fontsize);
legend_strs={'temperature','sum of squared error'}; 
if numel(vars_show)==1
ylabel(legend_strs(2),'FontSize',label_fontsize); 
else
    legend(legend_strs(vars_show));
end
grid on;
xlabel('number of iterations','FontSize',label_fontsize); % title('Parameter fitting by simulated annealing','FontSize',22)
% 
fig_subpl2=subplot(1,2,2); 
barplot_gca=barh(data_init_optim(:,sel_nodes)'); set(fig_subpl2,'ytick',1:numel(sel_nodes)); 
legend({'initial guess','data','optimized'},'FontSize',label_fontsize,'Box', 'off')
set(gca,'FontSize',label_ticks_fontsize); 
% xticklabels=get(gca,'xtick'); set(fig_subpl2,'xticklabel',xticklabels,'FontSize',label_ticks_fontsize);
xlabel('stationary probabilities','FontSize',label_fontsize); 
% set(fig_subpl2,'yticklabel',''); % strrep(nodes(sel_nodes),'_','\_'),'FontSize',label_fontsize
set(fig_subpl2,'yticklabel',strrep(nodes(sel_nodes),'_','\_'),'FontSize',label_fontsize);
