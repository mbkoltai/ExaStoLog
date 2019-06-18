function [stationary_state_vals_scan,stationary_node_vals_scan,stationary_state_inds_scan]=fcn_onedim_parscan_calc(stg_table,transition_rates_table,x0,nodes,...
                                                                                                parscan_matrix,scan_params,scan_params_up_down)
                                                                                            
[~,scan_par_inds,~]=fcn_get_trans_rates_tbl_inds(scan_params,scan_params_up_down,nodes);

if size(parscan_matrix,2)~=numel(scan_par_inds)
	error('parscan_matrix and scan_par_inds don''t have the same dimension, check them again!')
end

% stg_table_flag='stg_table_exists';
% if isempty(stg_table_flag)
%     [stg_table,~,~]=fcn_build_stg_table(truth_table_filename,nodes,transition_rates_table,'');
% end

stationary_node_vals_scan=zeros(numel(scan_par_inds),size(parscan_matrix,1),numel(nodes)); % transition_rates_table_mod=transition_rates_table;

scan_size = size(parscan_matrix,2);
for k=1:scan_size
    disp(strcat(num2str(round(100*k/scan_size)),'% done'))
    for val_counter=1:numel(parscan_matrix(:,k))
        transition_rates_table_mod=transition_rates_table; 
        transition_rates_table_mod(scan_par_inds(k))=parscan_matrix(val_counter,k);
        [A_sparse,~]=fcn_build_trans_matr(stg_table,transition_rates_table_mod,nodes,'');
        [stat_sol,~,~]=split_calc_inverse(A_sparse,transition_rates_table,x0);
        % sols per node
        [stationary_node_vals,~]=fcn_calc_init_stat_nodevals(x0,stat_sol);
        stationary_node_vals_scan(k,val_counter,:)=stationary_node_vals;
        if k==1 && val_counter==1
            stationary_state_vals_scan=zeros( numel(scan_par_inds),size(parscan_matrix,1), sum(stat_sol>0) );
            stationary_state_vals_scan(k,val_counter,:)=stat_sol(stat_sol>0);
            % index of nonzero states
            stationary_state_inds_scan = cell(numel(scan_par_inds),size(parscan_matrix,1));
            stationary_state_inds_scan{k,val_counter}=find(stat_sol>0);
        else
            stationary_state_vals_scan(k,val_counter,:)=stat_sol(stat_sol>0);
            stationary_state_inds_scan{k,val_counter}=find(stat_sol>0);
        end
    end
end

