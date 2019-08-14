% script for parameter fitting

% go to the folder of the file
% editor_service = com.mathworks.mlservices.MLEditorServices; editor_app = editor_service.getEditorApplication;
% active_editor = editor_app.getActiveEditor; storage_location = active_editor.getStorageLocation;
% file = char(storage_location.getFile); path_to_toolbox = fileparts(file); cd(path_to_toolbox);

% ADD FUNCTIONS to PATH
add_functions; model_name='breast_cancer_zanudo2017';
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
n_nodes=numel(nodes); truth_table_inputs=rem(floor([0:((2^n_nodes)-1)].'*pow2(0:-1:-n_nodes+1)),2);
x0=ones(2^n_nodes,1)/(2^n_nodes);
% get subgraphs
tic; [stat_sol,term_verts_cell,cell_subgraphs]=split_calc_inverse(A_sparse,transition_rates_table,x0); toc
% define x0 for given subgraph
output_nodes_zero=sum(truth_table_inputs(:,ismember(nodes,{'Apoptosis','Proliferation'})),2)==0;
subgr_ind=6;
pos_inds=ismember(1:numel(x0),cell_subgraphs{subgr_ind}) & output_nodes_zero'; % sum(pos_inds)
x0=zeros([2^numel(nodes),1]); x0(pos_inds)=1/sum(pos_inds); % sum(truth_table_inputs(x0>0,[1:3 7 8 18]))/sum(x0>0)

%% FITTING

writetable(table([find(stat_sol>0) stat_sol(stat_sol>0)]),'doc/sample_plots/breast_cancer_zanudo2017/zanudo2017_stat_sol.csv')

scan_params_filtered =[6 16 17]; scan_params_up_down_filtered = {1,1,[1 2]};
[~,~,predictor_names]=fcn_get_trans_rates_tbl_inds(scan_params_filtered,scan_params_up_down_filtered,nodes); 
% define data vector (generate some data OR load from elsewhere)
sel_param_vals=lognrnd(1,1,1,numel(predictor_names)); % abs(normrnd(1,0.5,1,numel(predictor_names)));
transition_rates_table=fcn_trans_rates_table(nodes,'uniform',[],[],predictor_names,sel_param_vals);
y_data=fcn_calc_init_stat_nodevals(x0,split_calc_inverse(fcn_build_trans_matr(stg_table,transition_rates_table,''),transition_rates_table,x0),'x0');

[fcn_statsol_sum_sq_dev,fcn_statsol_values]=fcn_simul_anneal(y_data,x0,stg_table,nodes,predictor_names);

% FITTING by simulated annealing (look at arguments in anneal/anneal.m)
% initial guess for parameters
init_vals=[0.01 10 10 0.01]; % lognrnd(0,2,size(predictor_names)); init_error=fcn_statsol_sum_sq_dev(init_vals);
% initial value of model nodes
y_init=fcn_calc_init_stat_nodevals(x0,...
    split_calc_inverse(fcn_build_trans_matr(stg_table,fcn_trans_rates_table(nodes,'uniform',[],[],predictor_names,init_vals),''),...
    transition_rates_table,x0),'');

fitting_arguments=struct('Verbosity',2, 'StopVal', 0.01);
[optim_par_vals,best_error,T_loss]=anneal(fcn_statsol_sum_sq_dev,init_vals,fitting_arguments);

% output with fitted parameters
y_optim_param=fcn_calc_init_stat_nodevals(x0,...
    split_calc_inverse(fcn_build_trans_matr(stg_table,fcn_trans_rates_table(nodes,'uniform',[],[],predictor_names,optim_par_vals),''),...
    transition_rates_table,x0),'');
% plot data, initial values, optimized values
data_init_optim=[y_data; y_init; y_optim_param]; min_val=min(min(data_init_optim(:,3:end))); max_val=max(max(data_init_optim(:,3:end)));

writetable(table(T_loss),'doc/sample_plots/breast_cancer_zanudo2017/zanudo2017_simul_ann_T_loss.csv')
writetable(table(data_init_optim),'doc/sample_plots/breast_cancer_zanudo2017/zanudo2017_data_init_optim.csv')