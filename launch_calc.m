% build STG
% tic; stg_table=fcn_build_stg_table(truth_table_filename,nodes); toc

% build A sparse
% tic; [A_sparse,~]=fcn_build_trans_matr(stg_table,transition_rates_table,''); toc

% generate initial conditions
% tic; x0=fcn_define_initial_states(initial_fixed_nodes,initial_fixed_nodes_vals,dom_prob,nodes,distrib_types{1},plot_flag); toc

% tic; [stat_sol,term_verts_cell,cell_subgraphs]=split_calc_inverse(A_sparse,transition_rates_table,x0); toc

%% looking inside 'split_calc_inverse' (subfunctions)

% crashes with models with large cyclic attractor(s)
% [vert_topol_sort,term_cycles_ind,~,~,term_cycle_bounds]=fcn_metagraph_scc(A_sparse_sub);

%% timing onedim parscan

% tic;
% [stationary_state_vals_onedimscan,stationary_node_vals_onedimscan,stationary_state_inds_scan]=...
%     fcn_onedim_parscan_calc(stg_table,transition_rates_table,x0,nodes,parscan_matrix,scan_params,scan_params_up_down);
% toc;

% 3 sec per parameter set, mostly consumed by regen of A_sparse
% after improvements down to 2.2 sec per param set with 20 node model

%% multidim parscan

% lhs_scan_dim=50; par_min_mean=-2; max_stdev=2; 
% [all_par_vals_lhs,stat_sol_nodes_lhs_parscan,stat_sol_states_lhs_parscan,stat_sol_states_lhs_parscan_cell]=... % outputs
%     fcn_multidim_parscan_latinhypcube(par_min_mean,max_stdev,sampling_type,lhs_scan_dim, ...
%                                             scan_params_sensit,scan_params_up_down_sensit, ... % transition rates
%                                             transition_rates_table,stg_table,x0,nodes);
          
% par_ind_table=[repelem(scan_params_sensit, cellfun(@(x) numel(x),scan_params_up_down_sensit))', horzcat(scan_params_up_down_sensit{:})'];
% trans_rate_scan_inds=(par_ind_table(:,1)-1)*2 + par_ind_table(:,2);
% %
% A_dim=2^numel(nodes);
% % states corresponding to the transition rates
% trans_matr_inds=cell2mat(arrayfun(@(x) find(ismember(2*(stg_table(:,3)-1)+stg_table(:,4),x)),trans_rate_scan_inds,'un',0)); % 0.8 sec
% trans_matr_inds_length=cell2mat(arrayfun(@(x) sum(ismember(2*(stg_table(:,3)-1)+stg_table(:,4),x)),trans_rate_scan_inds,'un',0)); % 0.8 sec
% n_par=numel(trans_rate_normalized);
% % IN THE LOOP
% tic;
% transition_rates_table_mod(trans_rate_scan_inds) = all_par_vals_lhs(k,:);
% trans_rate_normalized=transition_rates_table_mod(trans_rate_scan_inds)/sum(transition_rates_table_mod(:));
% norm_factor = sum(transition_rates_table_mod(:))/sum(transition_rates_table(:));
% A_sparse_mod=A_sparse/norm_factor; 
% % diagonal to 0, it'll be recalculated
% A_sparse_mod(1:(A_dim+1):numel(A_sparse_mod))=0;
% i_row=stg_table(trans_matr_inds,1); j_col=stg_table(trans_matr_inds,2); 
% seq_inds=i_row+(j_col-1)*A_dim;
% % reassign relevant trans rates
% A_sparse_mod(seq_inds)=cell2mat(arrayfun(@(x) repmat(trans_rate_normalized(x),trans_matr_inds_length(x),1), 1:n_par,'un',0)'); 
% % diagonal has to be recalculated
% A_sparse_mod = A_sparse_mod + speye(size(A_sparse_mod)) - diag(sum(A_sparse_mod,2)); % 0.33sec
% toc;
% % tic; [B,~]=fcn_build_trans_matr(stg_table,transition_rates_table_mod,''); toc % 1.7 sec


%% Sobol param scan

sobol_sensit_index=fcn_multidim_parscan_sobol_sensit_index([],var_types{1},...
                      all_par_vals_lhs,stat_sol_nodes_lhs_parscan,stat_sol_states_lhs_parscan,...
                      sample_size,... % # of calculations per parameter
                      sequential_indices_lhs,... % this is indices of transition rates in the original LHS
                      scan_params_filtered,scan_params_up_down_filtered,...% scan_params_sensit,scan_params_up_down_sensit
                      stg_table,transition_rates_table,x0,nodes,sel_nodes,plot_settings);

