function stat_sol_lhs_parscan=fcn_calc_paramsample_table(all_par_vals_lhs,scan_params,scan_params_up_down,transition_rates_table,stg_table,x0,nodes)


par_ind_table=[repelem(scan_params, cellfun(@(x) numel(x),scan_params_up_down))', horzcat(scan_params_up_down{:})'];
trans_rate_scan_inds=sub2ind(size(transition_rates_table),par_ind_table(:,2),par_ind_table(:,1));

transition_rates_table_mod=transition_rates_table;
stat_sol_lhs_parscan=zeros(size(all_par_vals_lhs,1),numel(nodes));
disp(strcat('dimension of parameter scan:',{' '},num2str(size(all_par_vals_lhs,1)),...
    {' '},'parameter sets of', {' '},num2str(size(all_par_vals_lhs,2)),{' '},'parameters.'))

lhs_scan_dim=size(all_par_vals_lhs,1);

for k=1:lhs_scan_dim

transition_rates_table_mod(trans_rate_scan_inds) = all_par_vals_lhs(k,:);
[A_sparse,~]=fcn_build_trans_matr(stg_table,transition_rates_table_mod,nodes,'');
[stat_sol,~,~]=split_calc_inverse(A_sparse,transition_rates_table_mod,x0);
[stationary_node_vals,~]=fcn_calc_init_stat_nodevals(x0,stat_sol);
stat_sol_lhs_parscan(k,:) = stationary_node_vals;

if rem(100*k/lhs_scan_dim,1)==0
    disp(strcat(num2str(round(100*k/lhs_scan_dim)),'% done'))
end

end
