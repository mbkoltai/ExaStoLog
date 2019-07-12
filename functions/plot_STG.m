%% STG plot of a model
function full_stg_plot=plot_STG(A_sparse,subgraph_index,term_verts_cell,cell_subgraphs,stat_sol,plot_settings,title_string,source_color)

if isempty(subgraph_index)
    A_digraph = digraph(A_sparse,'omitselfloops'); 
    term_verts_cell_subgraph = cell2mat(term_verts_cell);
    A_sparse_subgraph = A_sparse;
else
    A_digraph = digraph(A_sparse(cell_subgraphs{subgraph_index},cell_subgraphs{subgraph_index}),'omitselfloops'); 
    term_verts_cell_orig=cell2mat(term_verts_cell{subgraph_index});
    term_verts_cell_subgraph = find(ismember(cell_subgraphs{subgraph_index},cell2mat(term_verts_cell{subgraph_index})));
    A_sparse_subgraph=A_sparse(cell_subgraphs{subgraph_index},cell_subgraphs{subgraph_index});
end

max_term_verts_cell_subgraph = max(stat_sol(term_verts_cell_orig));

default_settings_cell=num2cell(plot_settings);
[fontsize,linewidth_val,arrowsize,default_markersize, source_sink_markersize]=deal(default_settings_cell{:});

full_stg_plot=plot(A_digraph,'Layout','force','ArrowSize',arrowsize,'MarkerSize',default_markersize,'LineWidth',linewidth_val);
xlim([min(full_stg_plot.XData) max(full_stg_plot.XData)]); ylim([min(full_stg_plot.YData) max(full_stg_plot.YData)]); 
% terminal nodes
for k=1:numel(term_verts_cell_subgraph)
    markersize_val=default_markersize+(source_sink_markersize-default_markersize)*stat_sol(term_verts_cell_orig(k))/max_term_verts_cell_subgraph;
    highlight(full_stg_plot,term_verts_cell_subgraph(k),'MarkerSize',markersize_val, 'NodeColor','red')
    hold off;
end

% source vertices
n_precision=3; 
source_vertices_subgraph=find(round(sum(A_sparse_subgraph - diag(diag(A_sparse_subgraph))),n_precision)==0);
highlight(full_stg_plot,source_vertices_subgraph,'MarkerSize',source_sink_markersize/2,'NodeColor',source_color); 


% label_strs = strsplit(num2str(find(diag(A)==1)')); labelnode(full_stg_plot,find(diag(A)==1)',label_strs)
if ~isempty(title_string)
    title(title_string, 'FontWeight','normal', 'FontSize', fontsize)
end
