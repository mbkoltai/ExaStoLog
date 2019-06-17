function [predictor_names,predictorImportance_vals]=fcn_multidim_parscan_predictorimport(scan_params,scan_params_up_down,all_par_vals_lhs,scan_values,...
                                                                                                nodes,sel_nodes,plot_type_flag)

par_ind_table=[repelem(scan_params, cellfun(@(x) numel(x),scan_params_up_down))', horzcat(scan_params_up_down{:})'];
prefix={'u_','d_'}; predictor_names=arrayfun(@(x) strcat(prefix(par_ind_table(x,2)), nodes(par_ind_table(x,1))), 1:size(par_ind_table,1),'un',0);
predictor_names=vertcat(predictor_names{:})';
% predictor importance values
predictorImportance_vals=cell2mat(arrayfun(@(x) ...
    predictorImportance(fitrtree(all_par_vals_lhs,scan_values(:,x),'PredictorNames',predictor_names)), sel_nodes,'un',0)');
max_y_val=max(predictorImportance_vals(:));
fontsize=14;

if strfind(plot_type_flag,'line')
% lineplot
if sum(sum(predictorImportance_vals==0))==0
    semilogy(1:numel(predictor_names), predictorImportance_vals,'LineWidth',2,'Marker','o'); 
else
    plot(1:numel(predictor_names), predictorImportance_vals,'LineWidth',2,'Marker','o'); 
end
set(gca,'XTick',1:numel(predictor_names)); set(gca,'xticklabel',strrep(predictor_names,'_','\_')); 
if size(scan_values,2)==numel(nodes)
    legend(nodes(sel_nodes),'FontWeight','normal','FontSize',fontsize,'Interpreter','none','Interpreter','none'); 
else
    legend(arrayfun(@(x) strcat('state_', num2str(sel_nodes(x))), 1:numel(sel_nodes),'un',0),'FontWeight','normal','FontSize',fontsize,'Interpreter','none')
end

elseif strfind(plot_type_flag,'bar')
% barplot
n_row_plot=round(sqrt(numel(sel_nodes))); n_col_plot=n_row_plot; if n_row_plot*n_col_plot<numel(sel_nodes); n_col_plot=n_row_plot+1; end
for k=1:numel(sel_nodes)
subplot(n_row_plot,n_col_plot,k); 
bar(predictorImportance_vals(k,:)); xlim([0 size(predictorImportance_vals,2)]+0.5)
if size(scan_values,2)==numel(nodes)
    title(nodes(sel_nodes(k)),'FontWeight','normal','FontSize',fontsize,'Interpreter','none'); 
else
%     n_nodes=numel(nodes); truth_table_inputs=rem(floor([0:((2^n_nodes)-1)].'*pow2(0:-1:-n_nodes+1)),2);
%     state_names=param_settings(3:end); suptitle(strcat('state #',num2str(state_names(var_ind)),', ',...
%         'p([',num2str( truth_table_inputs(state_names(var_ind),:) ),'])') );
title(strcat('state ',num2str(sel_nodes(k))),'FontWeight','normal','FontSize',fontsize,'Interpreter','none'); 
end


ylim([0 1.1*max_y_val]);
h=gca; 
if k>=numel(sel_nodes)-n_col_plot+1; h.XTickLabel=predictor_names; h.XTickLabelRotation = 45; h.TickLabelInterpreter = 'none';end 
% xlabel('Predictors'); 
if rem(k,n_col_plot)==1; ylabel('predictor importance'); end

end

else
    errorbar('Enter ''lineplot'' or ''barplot'' for <plot_type_flag> argument')    
end