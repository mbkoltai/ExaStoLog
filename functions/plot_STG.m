%% STG plot of a model
function full_stg_plot=plot_STG(A,counter,default_settings,xlim_vals,ylim_vals,titles,source_color)

% counter=4; fontsize=20;
A_digraph = digraph(A,'omitselfloops'); subnetws = conncomp(A_digraph,'Type','weak'); 

if length(unique(subnetws)) > 1 && ~isempty(counter)
    n_plot=2;
    subplot(1,n_plot,1); 
else
    n_plot=1;
end

default_settings_cell=num2cell(default_settings);
[fontsize,linewidth_val,arrowsize,default_markersize, source_sink_markersize]=deal(default_settings_cell{:});

full_stg_plot=plot(A_digraph,'Layout','force', 'ArrowSize', arrowsize, 'MarkerSize', default_markersize, 'LineWidth', linewidth_val); 
n_precision=3; terminal_nodes=find(round(diag(A),n_precision)==1);
highlight(full_stg_plot,terminal_nodes,'MarkerSize',source_sink_markersize, 'NodeColor','red')
source_vertices=find(round(sum(A - diag(diag(A))),n_precision)==0);
highlight(full_stg_plot,source_vertices,'MarkerSize',source_sink_markersize, 'NodeColor',source_color); 

if ~isempty(xlim_vals) && ~isempty(ylim_vals)
    xlim(xlim_vals(1,:)); ylim(ylim_vals(1,:)); 
end

% label_strs = strsplit(num2str(find(diag(A)==1)')); labelnode(full_stg_plot,find(diag(A)==1)',label_strs)
if ~isempty(titles)
    title(titles{1}, 'FontWeight','normal', 'FontSize', fontsize)
end

if n_plot>1 && ~isempty(counter)

A_i=A(subnetws==counter,subnetws==counter); % K_i=K(subnetws==counter,subnetws==counter);
subplot(1,n_plot,2); h=plot(digraph(A_i,'omitselfloops'),'Layout','force', 'ArrowSize', arrowsize); 
if ~isempty(xlim_vals) && ~isempty(ylim_vals)
    xlim(xlim_vals(2,:)); ylim(ylim_vals(2,:)); 
end
% xlim([-5 5]); ylim([-5 5]); 
% label_strs = strsplit(num2str(find(diag(A_i)==1)')); labelnode(h,find(diag(A_i)==1)',label_strs)
% highlight(h,1:size(A_i,1),'MarkerSize',2);
terminal_nodes = find(round(diag(A_i),n_precision)==1);
highlight(h,terminal_nodes,'MarkerSize',source_sink_markersize, 'NodeColor','red'); 
source_vertices=find(round(sum(A_i - diag(diag(A_i))),n_precision)==0);
highlight(h,source_vertices,'MarkerSize',source_sink_markersize, 'NodeColor','green'); 
title(titles{2},'FontWeight','normal','FontSize',fontsize)

end
