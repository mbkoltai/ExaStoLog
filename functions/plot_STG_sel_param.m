function []=plot_STG_sel_param(A_sparse,subgraph_index,term_verts_cell,cell_subgraphs,stat_sol,nodes,sel_params,sel_params_up_down,...
                                 stg_table,plot_settings,highlight_settings,tight_subpl,tight_subplot_pars)

% identify all parameters that play a role
sel_params_rep=cell2mat(arrayfun(@(x) repmat(sel_params(x),1,numel(sel_params_up_down{x})), 1:numel(sel_params_up_down),'un', 0));

sel_params_up_down = cell2mat(sel_params_up_down);
if size(sel_params_up_down,1)>1
    sel_params_up_down=sel_params_up_down';
end

selected_pars=sub2ind([numel(nodes) 2], sel_params_rep,sel_params_up_down);
trans_rate_indices=sub2ind([numel(nodes) 2], stg_table(:,3),stg_table(:,4));

if strcmp(selected_pars,'all')   
    relev_params=unique(trans_rate_indices)';
else
    relev_params=selected_pars;
end

rate_names={strcat('u_',nodes) strcat('d_',nodes)}; rate_names=horzcat(rate_names{:}); rate_names_sel=rate_names(relev_params);

if ~isempty(subgraph_index)
    if isempty(cell_subgraphs)
       disp('provide subgraph indices!');
    else
        rel_edges=cell_subgraphs{subgraph_index}; % A_sparse=A_sparse(cell_subgraphs{counter},cell_subgraphs{counter});
    end
else
    rel_edges=1:size(A_sparse,1);
end

n=numel(relev_params); round_sqrt=round(sqrt(n));
if round_sqrt^2 >= n
    n_row_plot=round_sqrt; n_col_plot=round_sqrt;
else
    n_row_plot=round_sqrt; n_col_plot=round_sqrt+1;
end

if ~isempty(tight_subpl)
    [ha,~] = tight_subplot(n_row_plot,n_col_plot,tight_subplot_pars(1,:),tight_subplot_pars(2,:),tight_subplot_pars(3,:));
end
% tight_subplot(Nh, Nw, gap, marg_h, marg_w): 
% gap: gaps between the axes in normalized units
% marg_h: margins in height in normalized units
% marg_w: margins in width in normalized units

for k=1:numel(relev_params)
    % subplot(ceil(sqrt(numel(relev_params))), ceil(sqrt(numel(relev_params))), k); 
    if isempty(tight_subpl)
        subplot(n_row_plot,n_col_plot,k);
    else
        axes(ha(k));
    end
    
    % rows and indices of entire A_sparse
    trans_matr_rows=stg_table(trans_rate_indices==relev_params(k),1); trans_matr_cols=stg_table(trans_rate_indices==relev_params(k),2); 
    % rows and cols of subgraphs
%     if ~isempty(subgraph_index)
%         trans_matr_rows=find(ismember(cell_subgraphs{subgraph_index}, trans_matr_rows)); 
%         trans_matr_cols=find(ismember(cell_subgraphs{subgraph_index}, trans_matr_cols)); 
%     end
    row_inds=ismember(rel_edges,trans_matr_rows); col_inds=ismember(rel_edges,trans_matr_cols);
    
    if sum(row_inds)==0
        disp(strcat(rate_names_sel{k},' has no effect on any transition in this subgraph!'))
    else
        
    full_stg_plot=plot_STG(A_sparse,subgraph_index,term_verts_cell,cell_subgraphs,stat_sol,plot_settings,[]);
    % highlight
    A_selected_trans=zeros(size(A_sparse(cell_subgraphs{subgraph_index},cell_subgraphs{subgraph_index}))); 
    % select the highlighted edges
    A_selected_trans(row_inds,col_inds) = A_sparse(cell_subgraphs{subgraph_index}(row_inds),cell_subgraphs{subgraph_index}(col_inds));
    highlight(full_stg_plot,digraph(A_selected_trans,'omitselfloops'),'EdgeColor', highlight_settings{1},'LineWidth',highlight_settings{2})
    title(strrep(rate_names_sel{k},'_','\_'),'FontWeight','normal','FontSize',14)
    end
end

