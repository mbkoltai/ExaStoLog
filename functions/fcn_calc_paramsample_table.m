function [stat_sol_lhs_parscan,stat_sol_states_lhs_parscan]=fcn_calc_paramsample_table(paramsample_table,...
                                                                multiscan_pars,multiscan_pars_up_down,...
                                                                transition_rates_table,stg_table,x0,disp_var)

par_ind_table=[repelem(multiscan_pars, cellfun(@(x) numel(x),multiscan_pars_up_down))', horzcat(multiscan_pars_up_down{:})'];
trans_rate_scan_inds=(par_ind_table(:,1)-1)*2 + par_ind_table(:,2);
transition_rates_table_mod=transition_rates_table;

A_dim=2^size(transition_rates_table,2);
% states corresponding to the transition rates
% trans_matr_inds=cell2mat(arrayfun(@(x) find(ismember(2*(stg_table(:,3)-1)+stg_table(:,4),x)),trans_rate_scan_inds,'un',0)); % 0.8 sec
trans_matr_inds_cell=arrayfun(@(x) find(ismember(2*(stg_table(:,3)-1)+stg_table(:,4),x)),trans_rate_scan_inds,'un',0);
% trans_matr_inds_length=cell2mat(arrayfun(@(x) sum(ismember(2*(stg_table(:,3)-1)+stg_table(:,4),x)),trans_rate_scan_inds,'un',0)); % 0.8 sec
% i_row=cellfun(@(x) stg_table(x,1), trans_matr_inds_cell,'un',0); j_col=cellfun(@(x) stg_table(x,2), trans_matr_inds_cell,'un',0);
seq_inds=cellfun(@(x) stg_table(x,1) + (stg_table(x,2)-1)*A_dim, trans_matr_inds_cell,'un',0); % i_row+(j_col-1)*A_dim;

% n_par=numel(trans_rate_scan_inds);
[A_sparse,~]=fcn_build_trans_matr(stg_table,transition_rates_table,'');
stg_sorting_cell=fcn_scc_subgraphs(A_sparse,x0);

stat_sol_states_lhs_parscan=cell(size(paramsample_table,1),1); % zeros(size(all_par_vals_lhs,1),sum(stat_sol>0));
stat_sol_lhs_parscan=zeros(size(paramsample_table,1),size(transition_rates_table,2));
disp(strcat('dimension of parameter scan:',{' '},num2str(size(paramsample_table,1)),...
    {' '},'parameter sets of', {' '},num2str(size(paramsample_table,2)),{' '},'parameters.'))

lhs_scan_dim=size(paramsample_table,1);

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
[stat_sol,~,~]=split_calc_inverse(A_sparse_mod,stg_sorting_cell,transition_rates_table_mod,x0);
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