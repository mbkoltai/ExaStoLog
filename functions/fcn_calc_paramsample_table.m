function [stat_sol_lhs_parscan,stat_sol_states_lhs_parscan]=fcn_calc_paramsample_table(paramsample_table,...
                                                                multiscan_pars,multiscan_pars_up_down,...
                                                                transition_rates_table,stg_cell,x0,disp_var)

par_ind_table=[repelem(multiscan_pars, cellfun(@(x) numel(x),multiscan_pars_up_down))', horzcat(multiscan_pars_up_down{:})'];
trans_rate_scan_inds=(par_ind_table(:,1)-1)*2 + par_ind_table(:,2);
transition_rates_table_mod=transition_rates_table;
n_nodes=size(transition_rates_table,2);

A_dim=2^size(transition_rates_table,2);
% states corresponding to the transition rates
% trans_matr_inds_cell=arrayfun(@(x) find(ismember(2*(stg_table(:,3)-1)+stg_table(:,4),x)),trans_rate_scan_inds,'un',0);
% seq_inds=cellfun(@(x) stg_table(x,1) + (stg_table(x,2)-1)*A_dim, trans_matr_inds_cell,'un',0); % i_row+(j_col-1)*A_dim;
plus_minus_inds = par_ind_table(:,2)'; plus_minus_inds(plus_minus_inds==2)=-1;
row_inds=arrayfun(@(x) stg_cell{par_ind_table(x,2),par_ind_table(x,1)},1:size(par_ind_table,1),'un',0);
col_inds=arrayfun(@(x) stg_cell{par_ind_table(x,2),par_ind_table(x,1)}+plus_minus_inds(x)*2^(n_nodes-par_ind_table(x,1)),1:size(par_ind_table,1),'un',0);
seq_inds=arrayfun(@(x) row_inds{x} + (col_inds{x}-1)*A_dim, 1:numel(col_inds),'un',0); 

% i_row=stg_cell{scan_par_inds(k)}; 
% if rem(scan_par_inds(k),2)>0; j_col=i_row+2^(numel(nodes) - ceil(scan_par_inds(k)/2)); 
% else; j_col=i_row-2^(numel(nodes) - ceil(scan_par_inds(k)/2)); 
% end
% seq_inds=i_row+(j_col-1)*A_dim;

% n_par=numel(trans_rate_scan_inds);
% [A_sparse,~]=fcn_build_trans_matr(stg_table,transition_rates_table,'');
[A_sparse,~]=fcn_build_trans_matr_stgcell(stg_cell,transition_rates_table,'');
stg_sorting_cell=fcn_scc_subgraphs(A_sparse,x0);

stat_sol_states_lhs_parscan=cell(size(paramsample_table,1),1); % zeros(size(all_par_vals_lhs,1),sum(stat_sol>0));
stat_sol_lhs_parscan=zeros(size(paramsample_table,1),size(transition_rates_table,2));
disp(strcat('dimension of parameter scan:',{' '},num2str(size(paramsample_table,1)),...
    {' '},'parameter sets of', {' '},num2str(size(paramsample_table,2)),{' '},'parameters.'))

lhs_scan_dim=size(paramsample_table,1);
rows_with_zero=unique(paramsample_table>0,'rows'); rows_with_zero=rows_with_zero(sum(rows_with_zero,2)<size(paramsample_table,2),:);
stg_sorting_cell_zeros=cell(1,size(rows_with_zero,1));

for k=1:lhs_scan_dim

transition_rates_table_mod(trans_rate_scan_inds)=paramsample_table(k,:);
trans_rate_normalized=transition_rates_table_mod(trans_rate_scan_inds)/sum(transition_rates_table_mod(:));
norm_factor=sum(transition_rates_table_mod(:))/sum(transition_rates_table(:));
A_sparse_mod=A_sparse/norm_factor; 
% diagonal to 0, it'll be recalculated
A_sparse_mod(1:(A_dim+1):numel(A_sparse_mod))=0;

% reassign relevant trans rates
for cell_cntr=1:numel(seq_inds)
    A_sparse_mod(seq_inds{cell_cntr}) = trans_rate_normalized(cell_cntr);
end
% diagonal has to be recalculated
A_sparse_mod=A_sparse_mod + speye(size(A_sparse_mod)) - diag(sum(A_sparse_mod,2)); % 0.33sec

% calculate solution
if sum(paramsample_table(k,:)==0)>0
    [~,row_match]=ismember(paramsample_table(k,:)>0, rows_with_zero);
    which_zero_ind=row_match(1);
    if isempty(stg_sorting_cell_zeros{which_zero_ind}) 
        stg_sorting_cell_zeros{which_zero_ind}=fcn_scc_subgraphs(A_sparse_mod,x0);
    end
    [stat_sol,~,~]=split_calc_inverse(A_sparse_mod,stg_sorting_cell_zeros{which_zero_ind},transition_rates_table_mod,x0);
else
    [stat_sol,~,~]=split_calc_inverse(A_sparse_mod,stg_sorting_cell,transition_rates_table_mod,x0);
end
disp(k)

[stationary_node_vals,~]=fcn_calc_init_stat_nodevals(x0,stat_sol,'');
stat_sol_lhs_parscan(k,:) = stationary_node_vals;
stat_sol_states_lhs_parscan{k} = stat_sol(stat_sol>0);

if ~isempty(disp_var)
if rem(disp_var*k/lhs_scan_dim,1)==0
    disp(strcat(num2str(round(100*k/lhs_scan_dim)),'% done'))
end
end

end

if numel(unique(cell2mat(arrayfun(@(x) numel(stat_sol_states_lhs_parscan{x}), 1:numel(stat_sol_states_lhs_parscan),'un',0))))==1
    stat_sol_states_lhs_parscan=cell2mat(stat_sol_states_lhs_parscan')';
end