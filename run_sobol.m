% script for parameter fitting

% go to the folder of the file
% editor_service = com.mathworks.mlservices.MLEditorServices; editor_app = editor_service.getEditorApplication;
% active_editor = editor_app.getActiveEditor; storage_location = active_editor.getStorageLocation;
% file = char(storage_location.getFile); path_to_toolbox = fileparts(file); cd(path_to_toolbox);

% ADD FUNCTIONS to PATH
add_functions; model_name='EMT_cohen_ModNet';
% where to save figures
plot_save_folder=strcat('doc/sample_plots/',model_name,'/');
[nodes,rules]=fcn_bnet_readin(strcat('model_files/',model_name,'.bnet')); 
% if yes, we generate a function file, which will create the model
truth_table_filename='fcn_truthtable.m'; fcn_write_logicrules(nodes,rules,truth_table_filename)
% from the model we generate the STG table, that is independent of values of transition rates
% this takes 20-30 seconds for 20 node model
stg_table=fcn_build_stg_table(truth_table_filename,nodes);
% RATES
chosen_rates=[]; chosen_rates_vals=[];
distr_type={'uniform','random'}; % <uniform> assigns a value of 1 to all params. other option: <random>
meanval=[]; sd_val=[]; % if 'random' is chosen, the mean and standard dev of a normal distrib has to be defined 
transition_rates_table=fcn_trans_rates_table(nodes,distr_type{1},meanval,sd_val,chosen_rates,chosen_rates_vals);

% build transition matrix A with parameter values
tic; [A_sparse,~]=fcn_build_trans_matr(stg_table,transition_rates_table,''); toc
% set up x0
initial_fixed_nodes={'ECMicroenv','DNAdamage','Metastasis','Migration','Invasion','EMT','Apoptosis'}; 
initial_fixed_nodes_vals=[1 1 zeros(1,5)];
dom_prob=1;
% if <random> the probability is randomly distributed among states, if <uniform> uniformly
distrib_types={'random','uniform'};
% if plot_flag non-empty, we get a bar plot of initial values
plot_flag='';
% function assigns a probability of <dom_prob> to the states with the fixed nodes having the defined values
tic; x0=fcn_define_initial_states(initial_fixed_nodes,initial_fixed_nodes_vals,dom_prob,nodes,distrib_types{1},plot_flag); toc
%%%%
tic; stg_sorting_cell=fcn_scc_subgraphs(A_sparse,x0); toc
% calculated stationary solution
tic; [stat_sol,term_verts_cell,cell_subgraphs]=split_calc_inverse(A_sparse,stg_sorting_cell,transition_rates_table,x0); toc

%% LHS 
scan_params_sensit=[13 14 15 16 17 19]; scan_params_up_down_sensit={[1 2],[1 2],[1 2],2,2,2};

if isfile('doc/sample_plots/EMT_cohen_ModNet/ECMicroenv_DNAdamage_11/all_par_vals_lhs.csv') && ...
        isfile('doc/sample_plots/EMT_cohen_ModNet/ECMicroenv_DNAdamage_11/stat_sol_nodes_lhs_parscan.csv')
% ~exist('all_par_vals_lhs','var') && ~exist('stat_sol_nodes_lhs_parscan','var')

stat_sol_nodes_lhs_parscan=table2array(readtable('doc/sample_plots/EMT_cohen_ModNet/ECMicroenv_DNAdamage_11/stat_sol_nodes_lhs_parscan.csv'));
all_par_vals_lhs=table2array(readtable('doc/sample_plots/EMT_cohen_ModNet/ECMicroenv_DNAdamage_11/all_par_vals_lhs.csv'));
    
else
sampling_types={'lognorm','linear','logunif'};
sampling_type=sampling_types{3};
% % <lhs_scan_dim>: number of param sets
lhs_scan_dim=500; par_min_mean=-2; max_stdev=2; % repmat(0.5,1,numel(cell2mat(scan_params_up_down(:))));
[all_par_vals_lhs,stat_sol_nodes_lhs_parscan,stat_sol_states_lhs_parscan,stat_sol_states_lhs_parscan_cell]=... % outputs
    fcn_multidim_parscan_latinhypcube(par_min_mean,max_stdev,sampling_type,lhs_scan_dim, ...
                     scan_params_sensit,scan_params_up_down_sensit, ... % transition rates
                     transition_rates_table,stg_table,x0,nodes);

writetable(table(stat_sol_states_lhs_parscan),'doc/sample_plots/EMT_cohen_ModNet/ECMicroenv_DNAdamage_11/stat_sol_states_lhs_parscan.csv')
writetable(table(stat_sol_nodes_lhs_parscan),'doc/sample_plots/EMT_cohen_ModNet/ECMicroenv_DNAdamage_11/stat_sol_nodes_lhs_parscan.csv')
writetable(table(all_par_vals_lhs),'doc/sample_plots/EMT_cohen_ModNet/ECMicroenv_DNAdamage_11/all_par_vals_lhs.csv')

end

%% Sobol sensitivity

scan_params_filtered=scan_params_sensit; scan_params_up_down_filtered=scan_params_up_down_sensit;
[~,sequential_indices_lhs,~] = fcn_get_trans_rates_tbl_inds(scan_params_sensit,scan_params_up_down_sensit,nodes);

sel_nodes=[]; % 3:numel(nodes); % scan_params_sensit;
sample_size=200; % if left empty, the sample size is half of the original param scan <all_par_vals_lhs>
% PLOT SETTINGS: [fontsize_plot,fontsize_axes,fontsize_title, min_color(optional), max_color(opt), progress_calcul_every_x_% (opt)];
plot_settings=[20 30 30 NaN NaN 10];
var_types={'node','state'}; % analysis for states or nodes
% to calculate Sobol total sensitivity we need <sample_size*numel(scan_params_up_down)> evaluations of the model
sobol_sensit_index=fcn_multidim_parscan_sobol_sensit_index([],var_types{1},...
                      all_par_vals_lhs,stat_sol_nodes_lhs_parscan,[],...
                      sample_size,... % # of calculations per parameter
                      sequential_indices_lhs,... % this is indices of transition rates in the original LHS
                      scan_params_filtered,scan_params_up_down_filtered,...% scan_params_sensit,scan_params_up_down_sensit
                      stg_table,transition_rates_table,x0,nodes,sel_nodes,plot_settings);

                  
writetable(table(sobol_sensit_index),'doc/sample_plots/EMT_cohen_ModNet/ECMicroenv_DNAdamage_11/sobol_sensit_index.csv')