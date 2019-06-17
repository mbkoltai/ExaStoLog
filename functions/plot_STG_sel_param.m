function []=plot_STG_sel_param(A,counter,nodes,cell_subgraphs,all_selected,state_transitions_inds,default_settings,...
                                                highlight_settings,limits,tight_subpl,tight_subplot_pars)

if strcmp(all_selected,'all')
% identify all parameters that play a role
    trans_rate_indices=sub2ind([numel(nodes) 2], state_transitions_inds(:,3),state_transitions_inds(:,4));
    relev_params=unique(trans_rate_indices)';
else
    trans_rate_indices=sub2ind([numel(nodes) 2], state_transitions_inds(:,3),state_transitions_inds(:,4));
    relev_params=all_selected;
end

rate_names={strcat('u_',nodes) strcat('d_',nodes)}; rate_names=horzcat(rate_names{:}); rate_names_sel=rate_names(relev_params);

if ~isempty(counter)
    if isempty(cell_subgraphs)
       disp('provide subgraph indices!');
    else
        rel_edges=cell_subgraphs{counter}; A=A(cell_subgraphs{counter},cell_subgraphs{counter});
    end
else
    rel_edges=1:size(A,1);
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

if ~isempty(limits)
    xlim_vals=limits(1,:); ylim_vals=limits(2,:);
else
    xlim_vals=[]; ylim_vals=[]; 
end

for k=1:numel(relev_params)
    % subplot(ceil(sqrt(numel(relev_params))), ceil(sqrt(numel(relev_params))), k); 
    if isempty(tight_subpl)
        subplot(n_row_plot,n_col_plot,k);
    else
        axes(ha(k));
    end
    
    trans_matr_rows=state_transitions_inds(trans_rate_indices==relev_params(k),1); row_inds=ismember(rel_edges,trans_matr_rows);
    trans_matr_cols=state_transitions_inds(trans_rate_indices==relev_params(k),2); col_inds=ismember(rel_edges,trans_matr_cols);
    if sum(row_inds)==0
        disp(strcat(rate_names_sel{k},' has no effect on any transition in this subgraph!'))
    else
        
    full_stg_plot=plot_STG(A,'',default_settings,xlim_vals,ylim_vals,'','blue');
    % highlight
    A_selected_trans = zeros(size(A)); 
    % select the highlighted edges
    % disp(col_inds)
    A_selected_trans(row_inds,col_inds) = A(row_inds,col_inds);
    highlight(full_stg_plot,digraph(A_selected_trans,'omitselfloops'),'EdgeColor', highlight_settings{1},'LineWidth',highlight_settings{2}) 
    title(strrep(rate_names_sel{k},'_','\_'),'FontWeight','normal','FontSize',14)
    end
end

