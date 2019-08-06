
% build STG
% tic; stg_table=fcn_build_stg_table(truth_table_filename,nodes); toc

% build A sparse
% tic; [A_sparse,~]=fcn_build_trans_matr(stg_table,transition_rates_table,''); toc

% generate initial conditions
% tic; x0=fcn_define_initial_states(initial_fixed_nodes,initial_fixed_nodes_vals,dom_prob,nodes,distrib_types{1},plot_flag); toc

tic; [stat_sol,term_verts_cell,cell_subgraphs]=split_calc_inverse(A_sparse,transition_rates_table,x0); toc

