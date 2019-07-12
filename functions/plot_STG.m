%% STG plot of a model
function full_stg_plot=plot_STG(A_sparse,default_settings,xlim_vals,ylim_vals,titles,source_color)

A_digraph = digraph(A_sparse,'omitselfloops'); 

default_settings_cell=num2cell(default_settings);
[fontsize,linewidth_val,arrowsize,default_markersize, source_sink_markersize]=deal(default_settings_cell{:});

full_stg_plot=plot(A_digraph,'Layout','force', 'ArrowSize', arrowsize, 'MarkerSize', default_markersize, 'LineWidth', linewidth_val); 
n_precision=3; 
% terminal nodes
terminal_nodes=find(round(diag(A_sparse),n_precision)==1);
highlight(full_stg_plot,terminal_nodes,'MarkerSize',source_sink_markersize, 'NodeColor','red')
source_vertices=find(round(sum(A_sparse - diag(diag(A_sparse))),n_precision)==0);
highlight(full_stg_plot,source_vertices,'MarkerSize',source_sink_markersize, 'NodeColor',source_color); 

if ~isempty(xlim_vals) && ~isempty(ylim_vals)
    xlim(xlim_vals); ylim(ylim_vals); 
end

% label_strs = strsplit(num2str(find(diag(A)==1)')); labelnode(full_stg_plot,find(diag(A)==1)',label_strs)
if ~isempty(titles)
    title(titles{1}, 'FontWeight','normal', 'FontSize', fontsize)
end
