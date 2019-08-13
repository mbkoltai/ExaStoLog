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

% sobol_sensit_index=fcn_multidim_parscan_sobol_sensit_index([],var_types{1},...
%                       all_par_vals_lhs,stat_sol_nodes_lhs_parscan,stat_sol_states_lhs_parscan,...
%                       sample_size,... % # of calculations per parameter
%                       sequential_indices_lhs,... % this is indices of transition rates in the original LHS
%                       scan_params_filtered,scan_params_up_down_filtered,...% scan_params_sensit,scan_params_up_down_sensit
%                       stg_table,transition_rates_table,x0,nodes,sel_nodes,plot_settings);

%% open a figure

uiopen('doc/sample_plots/breast_cancer_zanudo2017/A_sol_init_101_PIMfree.fig',1)

%% init cond

node_inds=ismember(nodes,{'Alpelisib','Everolimus','PI3K','PIM','PDK1','PRAS40','TSC','mTORC1','Proliferation','Apoptosis'});

subgr_ind=find(ismember(input_nodes_per_subgraph, [0 1 0 0],'rows')); if numel(subgr_ind)>1; subgr_ind=subgr_ind(1); end
output_nodes_zero=sum(truth_table_inputs(:,ismember(nodes,{'Apoptosis','Proliferation'})),2)==0;
pos_inds=ismember(1:numel(x0),cell_subgraphs{subgr_ind}) & output_nodes_zero'; % sum(pos_inds)
x0=zeros(size(x0)); x0(pos_inds)=1/sum(pos_inds); % sum(truth_table_inputs(x0>0,[1:3 7 8 18]))/sum(x0>0)

% nodes(node_inds)
sel_nodes_inds=3:(numel(nodes)-2);
t=array2table(sum(truth_table_inputs(x0>0,sel_nodes_inds))/sum(x0>0)); t.Properties.VariableNames=nodes(sel_nodes_inds); t

%% STAT SOL

tic; [stat_sol,term_verts_cell,cell_subgraphs]=split_calc_inverse(A_sparse,transition_rates_table,x0); toc
% nonzero states displayed:
stat_sol(stat_sol>0)' % probability values of nonzero states
truth_table_inputs(stat_sol>0,:) % list of logical states that are nonzero
[stationary_node_vals,init_node_vals]=fcn_calc_init_stat_nodevals(x0,stat_sol,'x0');

sel_nodes_inds=3:numel(nodes);
t=array2table(truth_table_inputs(stat_sol>0,sel_nodes_inds)); t.Properties.VariableNames=nodes(sel_nodes_inds); t

%% PARCAN

parscan_min_max = [1e-2 1e2]; n_steps=2; sampling_types={'log','linear'}; 
parscan_matrix=fcn_onedim_parscan_generate_matrix(scan_params,scan_params_up_down,nodes,sampling_types{1},parscan_min_max,n_steps);

tic;
[stationary_state_vals_onedimscan,stationary_node_vals_onedimscan,stationary_state_inds_scan]=...
    fcn_onedim_parscan_calc(stg_table,transition_rates_table,x0,nodes,parscan_matrix,scan_params,scan_params_up_down);
toc;

%% PLOT

sensit_cutoff=0.01; 

figure()
[resp_coeff,scan_params_sensit,scan_params_up_down_sensit,fig_filename]=fcn_onedim_parscan_plot_parsensit(plot_types,plot_type_options,...
                                                   stationary_node_vals_onedimscan,stationary_state_vals_onedimscan,...
                                                   nonzero_states_inds,parscan_matrix,nodes,...
                                                   scan_params,scan_params_up_down,... % transition rates to scan in
                                                   sensit_cutoff,plot_param_settings);
