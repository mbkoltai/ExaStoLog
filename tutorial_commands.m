% commands for tutorial

% unzip the file toolboxes.zip!

% go to the folder of the file
editor_service=com.mathworks.mlservices.MLEditorServices; editor_app=editor_service.getEditorApplication;
active_editor=editor_app.getActiveEditor; storage_location=active_editor.getStorageLocation;
file=char(storage_location.getFile); path_to_toolbox=fileparts(file); cd(path_to_toolbox);

% ADD FUNCTIONS to PATH: required toolboxes are in <toolboxes.zip>
add_functions

%% READ IN model

% names of models
model_name_list = {'mammalian_cc', ...
'krasmodel15vars', ...
'breast_cancer_zanudo2017'....
'EMT_cohen_ModNet',...
'sahin_breast_cancer_refined'}; %
% select the index of one model
model_index=4;
model_name=model_name_list{model_index};

% read in model from BOOLNET file
[nodes,rules]=fcn_bnet_readin(strcat('model_files/',model_name,'.bnet'));
% Compare if rules and nodes are consistent
fcn_nodes_rules_cmp(nodes,rules)

% where to save figures
plot_save_folder=strcat('doc/sample_plots/',model_name,'/');
% write file with logical rules
truth_table_filename='fcn_truthtable.m'; fcn_write_logicrules(nodes,rules,truth_table_filename)
% build STG
tic; stg_table=fcn_build_stg_table(truth_table_filename,nodes); toc
% density of STG: size(stg_table,1)/(2^(2*numel(nodes)))

%% choose transition rates

chosen_rates=[]; chosen_rates_vals=[];
% ARGUMENTS
% <uniform> assigns a value of 1 to all params. <random> samples from a lognormal distribution
distr_type={'uniform','random'}; 
% if 'random' is chosen, the mean and standard dev of a normal distrib has to be defined
meanval=[]; sd_val=[]; 
transition_rates_table=fcn_trans_rates_table(nodes,distr_type{1},meanval,sd_val,chosen_rates,chosen_rates_vals);
                    
%% BUILD transition matrix

[A_sparse,~]=fcn_build_trans_matr(stg_table,transition_rates_table,'');

% visualize
% spy(A_sparse); xlabel('model states'); ylabel('model states'); set(gca,'FontSize',24)

%% Defining initial conditions

% selected nodes for inital conditions
initial_fixed_nodes_list=...
{ {'CycE','CycA','CycB','Cdh1','Rb_b1','Rb_b2','p27_b1','p27_b2'}, ... % mammalian_cc
{'cc','KRAS','DSB','cell_death'}, ... % krasmodel15vars
{'Alpelisib', 'Everolimus','PIM','Proliferation','Apoptosis'},...  % breast_cancer_zanudo2017
{'ECMicroenv','DNAdamage','Metastasis','Migration','Invasion','EMT','Apoptosis','Notch_pthw','p53'}, ... % EMT
{'EGF','ERBB1','ERBB2','ERBB3','p21','p27'}}; % sahin_breast_cancer_refined

% values for selected nodes
initial_fixed_nodes_vals_list = {[0 0 0 1 1 1 1 1], ... % mammalian_cc
            [1 1 1 0], ... % krasmodel15vars: [1 1] is cell cycle ON, KRAS mutation ON
            [0 1 0 zeros(1,2)],...  % breast_cancer_zanudo2017
            [1 1 zeros(1,5) 1 0],... % EMT-Cohen model: [0/1 0/1 zeros(1,5)]
            [1 0 0 0 1 1]}; % 1 zeros(1,numel(initial_fixed_nodes_list{model_index})-3) 1 1


% select the initial condition for the model we are working on
initial_fixed_nodes=initial_fixed_nodes_list{model_index}; 
initial_fixed_nodes_vals=initial_fixed_nodes_vals_list{model_index};

dom_prob=1;
% if <random> the probability is randomly distributed among states, if <uniform> uniformly
distrib_types={'random','uniform'};
% if plot_flag non-empty, we get a bar plot of initial values
plot_flag='';
% function assigns a probability of <dom_prob> to states with the fixed nodes
x0=fcn_define_initial_states(initial_fixed_nodes,initial_fixed_nodes_vals,...
				dom_prob,nodes,distrib_types{1},plot_flag);
            
%% Topological sorting of STG and identification of cycles

stg_sorting_cell=fcn_scc_subgraphs(A_sparse,x0);

%% calculate stationary solution

tic; [stat_sol,term_verts_cell,cell_subgraphs]=split_calc_inverse(A_sparse,stg_sorting_cell,transition_rates_table,x0); toc
% by model variables
[stationary_node_vals,init_node_vals]=fcn_calc_init_stat_nodevals(x0,stat_sol);

%% visualize solution

% ARGUMENTS
% matrix_input: [], K_sparse or A_sparse (kinetic or transition matrix)
% min_max_col: minimum and max color for heatmap
% fontsize: [font size on the heatmap, title font size for stationary solutions]
% barwidth_states_val: width of the bars for bar plot of stationary solutions of states
% sel_nodes: nodes to show. If left empty, all nodes are shown
% prob_thresh: minimal value for probability to display 
% (useful for visibility if many attractor states or large cyclic attractor(s))

% Call the function by:
sel_nodes=[];
min_max_col=[0 1]; barwidth_states_val=0.8;
% fontsize: [fontsize of plot, fontsize of titles, fontsize of binary states]
fontsize=[24 40 20]; 
plot_settings = [fontsize barwidth_states_val min_max_col]; prob_thresh=0.03;

figure('name','A_K_stat_sol')
fcn_plot_A_K_stat_sol(A_sparse,nodes,sel_nodes,stat_sol,x0,plot_settings,prob_thresh)

%% SAVE figure
if exist(plot_save_folder,'dir')==0; mkdir(plot_save_folder); end
fig_file_type={'.png','.eps','.pdf','.jpg','.tif'};
% if <overwrite_flag> non-empty then existing file with same name is overwritten.
overwrite_flag='yes';

% resolution of the figures (dpi)
resolution_dpi='-r350';
% SAVE
fcn_save_fig('single_solution_states_nodes_stat_sol',plot_save_folder,fig_file_type{3},overwrite_flag,resolution_dpi)

%% Visualize binary heatmap of nonzero stationary states

% stat_sol: vector of stationary solutions
% prob_thresh: probability threshold for states to show (if left empty, all states shown)
prob_thresh=0.01;  
% term_verts_cell: index of subgraphs for stable states
% nodes: name of nodes
% sel_nodes: nodes to show. if none selected, all nodes will be shown
sel_nodes=[];
% plot_param_settings
% num_size_plot: font size of 0/1s on the heatmap
% hor_gap: horizontal gap between terminal SCCs, bottom_marg: bottom margin, left_marg: left margin
numsize_plot=26; fontsize=36; hor_gap=0.02; bottom_marg=0.31; left_marg=0.22;
plot_param_settings=[numsize_plot fontsize hor_gap bottom_marg left_marg];
% tight_subplot_flag: want to use tight subplot? | ranking_flag: order states by probability?
tight_subplot_flag='yes'; ranking_flag='yes';

% PLOT
figure('name','statsol_binary_heatmap')
statsol_binary_heatmap=fcn_plot_statsol_bin_hmap(stat_sol,prob_thresh,...
term_verts_cell,nodes,sel_nodes,plot_param_settings,tight_subplot_flag,ranking_flag);
%% SAVE
resolution_dpi='-r350';
fcn_save_fig('binary_heatmap_states',plot_save_folder,fig_file_type{3},overwrite_flag,resolution_dpi);

%% select transition rates for 1-dimensional parameter scan

popul_subgraphs=cellfun(@(x) sum(ismember(find(x0>0), x)), cell_subgraphs)>0;
subgraph_states=cell2mat(cell_subgraphs(popul_subgraphs)');
par_inds_table=unique(stg_table(ismember(stg_table(:,1), subgraph_states) | ...
					ismember(stg_table(:,2), subgraph_states),3:4),'rows');
                
% most common transitions
for k=1:size(par_inds_table,1)
    param_freq(k) = sum(stg_table(:,3)==par_inds_table(k,1) & stg_table(:,4)==par_inds_table(k,2));
end
% top n most frequent transitions
[~,top_freq_trans_rates]=maxk(param_freq,6);

% all
scan_params=unique(par_inds_table(:,1))';
% most common: scan_params=par_inds_table(top_freq_trans_rates,1)';
% by name: scan_params=find(ismember(nodes,{'Notch_pthw','p53','EMTreg','FOXO3','p63_73'}));

% up and down rates (all)
scan_params_up_down=arrayfun(@(x) par_inds_table(par_inds_table(:,1)==x,2)', scan_params,'un',0);

%% CALCULATE 1-dimensional scan

% min and max of range of values; resolution of the scan; linear or logarithmic sampling
parscan_min_max = [1e-2 1e2]; n_steps=10; sampling_types={'log','linear'};

% matrix of parameter values
parscan_matrix=fcn_onedim_parscan_generate_matrix(scan_params,scan_params_up_down,...
		nodes,sampling_types{1},parscan_min_max,n_steps);

% calculation
[stationary_state_vals_onedimscan,stationary_node_vals_onedimscan,stationary_state_inds_scan]=...
    fcn_onedim_parscan_calc(stg_table,transition_rates_table,x0,nodes,parscan_matrix,scan_params,scan_params_up_down);

%% PLOT 1-dimensional scan grouped by transition rates

% index of nonzero states
nonzero_states_inds=find(stat_sol>0);
% plot parameters: [vertical, horizontal] gap between subplots, margins at [bottom, top], [left, right]
height_width_gap=[0.08 0.03]; bott_top_marg =[0.05 0.05]; left_right_marg=[0.04 0.01];
params_tight_subplots={height_width_gap bott_top_marg left_right_marg};
% plot_param_settings: [fontsize_axes,fs_title,fs_legend,linewidth,params_tight_subplots,model_name]
plot_param_settings={24,34,24,4,{height_width_gap bott_top_marg left_right_marg},model_name};
% plotting stater or variables (nodes)?
state_or_node_flags={'nodes','states'};
% cutoff for minimal variation to show a variable
diff_cutoff=0.15;
figure('name','onedim parscan by param')
[fig_filename,onedim_paramscan_output_cell]=fcn_onedim_parscan_plot_by_params(state_or_node_flags{1},...
		stationary_node_vals_onedimscan,stationary_state_vals_onedimscan,...
		nonzero_states_inds,parscan_matrix,nodes,...
		scan_params,scan_params_up_down,... % selected parameters
		diff_cutoff,... % minimal variation for variable to be shown on plot
		plot_param_settings);
    
%% PLOT 1-dimensional scan grouped by model variables/states

% nonzero states of the model
nonzero_states_inds=find(stat_sol>0);
% sensit_cutoff: minimal value for local sensitivity or variation of model/state values
sensit_cutoff=0.1; 
% parameters of plot
height_width_gap=[0.1 0.04]; bott_top_marg=[0.03 0.1]; left_right_marg=[0.07 0.02]; 
params_tight_subplots={height_width_gap bott_top_marg left_right_marg};
% plot_param_settings: [fontsize_axes,fontsize_title,params_tight_subplots,model_name]
plot_param_settings={30,30,params_tight_subplots,model_name,'colorbar'};
% plot_param_settings={12,14,[],model_name}; 
% select type of plot
plot_types={{'lineplot','heatmap'} {'nodes','states'} {'values','sensitivity'}};
plot_type_options=[1 2 1];

[resp_coeff,scan_params_sensit,scan_params_up_down_sensit,fig_filename]=...
	fcn_onedim_parscan_plot_parsensit(plot_types,plot_type_options,...
		stationary_node_vals_onedimscan,stationary_state_vals_onedimscan,...
		nonzero_states_inds,parscan_matrix,nodes,...
		scan_params,scan_params_up_down,...
		sensit_cutoff,plot_param_settings);
%% SAVE
fcn_save_fig(strcat(fig_filename,'_cutoff',strrep(num2str(sensit_cutoff),'.','p')),...
	plot_save_folder,fig_file_type{1},'overwrite','-r200');

%% PLOT local sensitivities

plot_type_options=[2 2 2];
figure('name',strjoin(arrayfun(@(x) plot_types{x}{plot_type_options(x)}, ...
		1:numel(plot_type_options), 'un',0),'_'));
[resp_coeff,scan_params_sensit,scan_params_up_down_sensit,fig_filename]=...
			fcn_onedim_parscan_plot_parsensit(plot_types,plot_type_options,...
                          stationary_node_vals_onedimscan,stationary_state_vals_onedimscan,...
                          nonzero_states_inds,parscan_matrix,nodes,...
                          scan_params,scan_params_up_down,... 
                          sensit_cutoff,plot_param_settings);
                      
%% SAVE
fcn_save_fig(strcat(fig_filename,'_cutoff',strrep(num2str(sensit_cutoff),'.','p')),...
	plot_save_folder,fig_file_type{1},'overwrite','-r200');

%% Multidimensional (2D) parameter scanning with regular grids

n_scanvals=10; scanvals=logspace(-2,2,n_scanvals);  % with zero: [0 logspace(-2,2,n_scanvals-1)]
meshgrid_scanvals=meshgrid(scanvals,scanvals);

paramsample_table=[repelem(scanvals,n_scanvals)' ...
	reshape(reshape(repelem(scanvals,n_scanvals),n_scanvals,n_scanvals)',n_scanvals^2,1)]; 
% transition rates to scan in
multiscan_pars=[11 13]; multiscan_pars_up_down={1 1};

disp_var=5; % show at every n% the progress
[stat_sol_paramsample_table,stat_sol_states_paramsample_table]=...
    fcn_calc_paramsample_table(paramsample_table,multiscan_pars,...
		multiscan_pars_up_down,transition_rates_table,stg_table,x0,disp_var);
    
%% PLOT as heatmap

% what model variables to plot?
sel_nodes=4;
% plot_settings: [fontsize on plot, fs axes, fs subplot titles, fs axes labels]
plot_settings=[28 30 40]; figure('name','2D scan')
fcn_plot_twodim_parscan(stat_sol_paramsample_table,scanvals,...
			multiscan_pars,multiscan_pars_up_down,...
			nodes,sel_nodes,plot_settings)

%% SAVE PLOT
resolution_dpi='-r200'; 
file_name_prefix=strcat('twodim_parscan_',strjoin(nodes(sel_nodes),'_'));
fcn_save_fig(file_name_prefix,plot_save_folder,fig_file_type{1},'overwrite',resolution_dpi);

%% Multidimensional parameter scanning with Latin Hypercube Sampling

sampling_types={'lognorm','linear','logunif'}; sampling_type=sampling_types{3};
% par_min_mean: minimum or (if lognormal) mean of distribution. 
%  				Can be a scalar or vector (if different values for different parameters)
% max_stdev: maximum or in case of lognormal the mean of distribution. Scalar or vector
%
% for 'lognorm','logunif' provide LOG10 value of desired mean/min & stdev/max (-2 is a mean of 0.01)
par_min_mean=-2; % repmat(1.5,1,numel(cell2mat(scan_params_up_down_sensit(:)))); par_min_mean(4)=3; 
max_stdev=2; 	 % repmat(0.5,1,numel(cell2mat(scan_params_up_down_sensit(:))));
% number of param sets
lhs_scan_dim=1000;

% RUN the LHS
[all_par_vals_lhs,stat_sol_nodes_lhs_parscan,stat_sol_states_lhs_parscan]=... % outputs
    fcn_multidim_parscan_latinhypcube(par_min_mean,max_stdev,sampling_type,lhs_scan_dim, ...
                          scan_params_sensit,scan_params_up_down_sensit, ...
                          transition_rates_table,stg_table,x0,nodes);
                      
%% Visualize multi-dimensional parameter scans by scatter plots

% which variable to plot?
var_ind=4;
% STATES or NODES? <scan_values>: values to be plotted
% model variables: stat_sol_nodes_lhs_parscan; states: stat_sol_states_lhs_parscan
scan_values=stat_sol_states_lhs_parscan; 

sampling_type=sampling_types{3}; % sampling_types={'lognorm','linear','logunif'};
% file_name_prefix=strcat('LHS_parscan_scatterplot_trend_',nodes{var_ind}); 
file_name_prefix=strcat('LHS_parscan_scatterplot_trend_state',num2str(var_ind));
% param_settings: [number_bins_for_mean,trendline_width,axes_fontsize,index nonzero states]
param_settings = [50 6 24 size(stat_sol_states_lhs_parscan)];

figure('name',num2str(var_ind))
fcn_multidim_parscan_scatterplot(var_ind,all_par_vals_lhs,scan_values,...
        scan_params_sensit,scan_params_up_down_sensit,nodes,sampling_type,param_settings)

%% SAVE
resolution_dpi='-r200'; 
fcn_save_fig(file_name_prefix,plot_save_folder,fig_file_type{3},'overwrite',resolution_dpi);

%% PLOT Correlations between model variables

% sel_nodes: name of selected nodes (pls provide in ascending order) (if left empty, all shown)
sel_nodes=[3 7 8 10 11 13:15 17:20]; 
% plot_settings: [fontsize on plot, fontsize on axes/labels]
plot_settings=[NaN 26 32]; 
% we'll plot correlations between variables, as a heatmap
plot_type_flag={'var_var','heatmap'}; 

figure('name',strjoin(plot_type_flag))
[varvar_corr_matr,p_matrix_vars]=...
	fcn_multidim_parscan_parvarcorrs(plot_type_flag,all_par_vals_lhs,...
		stat_sol_nodes_lhs_parscan,nodes,sel_nodes,[],[],[],plot_settings);

%% SAVE
resolution_dpi='-r350'; 
fcn_save_fig(fig_prefix,plot_save_folder,fig_file_type{3},'overwrite',resolution_dpi);

%% Regression of variables by transition rates

plot_type_flag={'par_var','heatmap','r_sq'}; % {'par_var','heatmap'/'lineplot','r_sq'/'slope'}
sel_nodes=[];
% plot_settings=[fontsize,maximum value for heatmap colors], 
% if plot_settings(3)=NaN, then max color automatically selected
plot_settings=[30 30 0.29]; 
% if regression type is 'linlog', then the fit is y = a + b*log10(x)
regr_types={'log','linear'}; % log recommended if parameter values log-uniformly distributed in sampling
figure('name',strjoin(plot_type_flag))
scan_values=stat_sol_states_lhs_parscan; % or: stat_sol_nodes_lhs_parscan
[r_squared,slope_intercept]=fcn_multidim_parscan_parvarcorrs(plot_type_flag,...
		all_par_vals_lhs,scan_values,...
		nodes,sel_nodes,... % which nodes
		scan_params_sensit,scan_params_up_down_sensit, ... % same params as in LHS!
		regr_types{1},plot_settings);

%% SAVE plot
fig_prefix=strjoin(plot_type_flag,'_'); resolution_dpi='-r350'; 
fcn_save_fig(fig_prefix,plot_save_folder,fig_file_type{3},'overwrite',resolution_dpi)

%% Sobol total sensitivity index

% for indexing, we need the sequential indices of transition rates 
% (eg. 5th node's up rate {5,1}->9, 6th node down rate is {6,2}->12)
[par_ind_table,sequential_indices_lhs,~] = ...
	fcn_get_trans_rates_tbl_inds(scan_params_sensit,scan_params_up_down_sensit,nodes);
% threshold for R^2
r_sq_thresh=0.05;
% select transition rates
par_ind_table_filtered=par_ind_table(sum(r_squared>r_sq_thresh)>0,:);
scan_params_filtered=unique(par_ind_table_filtered(:,1))'; 
scan_params_up_down_filtered=...
	arrayfun(@(x) par_ind_table_filtered(par_ind_table_filtered(:,1)==x,2)',scan_params_filtered,'un',0);

% CALCULATION
sel_nodes=[]; 	% if left empty, all nodes/states are analyzed
sample_size=[]; % if left empty, the sample size is half of the original param scan <all_par_vals_lhs>
% how often (what %) should the progress of calculation be displayed?
disp_freq=10;
var_types={'node','state'}; % analysis for states or nodes
sobol_sensit_index=fcn_multidim_parscan_sobol_sensit_index([],var_types{2},...
	all_par_vals_lhs,stat_sol_nodes_lhs_parscan,stat_sol_states_lhs_parscan,...
	sample_size,... % # of calculations per parameter
	sequential_indices_lhs,... % indices of transition rates in the original LHS
	scan_params_filtered,scan_params_up_down_filtered,... % or: scan_params_sensit,scan_params_up_down_sensit
	stg_table,transition_rates_table,x0,nodes,sel_nodes,plot_settings,disp_freq);

%% PLOT

% plot_settings: 
% [fontsize_plot,fs_axes,fs_title,min_color(optional),max_color(opt),angle x-axis labels];
plot_settings=[30 30 40 0 0.5 90];
fcn_multidim_parscan_sobol_sensit_index(sobol_sensit_index,var_types{2},all_par_vals_lhs,[],[],[],...
	sequential_indices_lhs,scan_params_filtered,scan_params_up_down_filtered,[],[],[],...
	nodes,sel_nodes,plot_settings,[]);
xticklabels({'Metastasis','Apoptosis (p53)','Apoptosis (p63_73)'})

%% SAVE
resolution_dpi='-r350'; 
fcn_save_fig('sobol_sensitivity_index',plot_save_folder,fig_file_type{3},'overwrite',resolution_dpi)

%% Simulated annealing

% names of selected transition rates
[~,~,predictor_names]=fcn_get_trans_rates_tbl_inds(scan_params_sensit,scan_params_up_down_sensit,nodes); 

% define data vector (generate some data OR load from elsewhere)
data_param_vals=lognrnd(1,1,1,numel(predictor_names)); 
transition_rates_table_optim=fcn_trans_rates_table(nodes,'uniform',[],[],predictor_names,data_param_vals);

y_data=fcn_calc_init_stat_nodevals(x0,...
	split_calc_inverse(fcn_build_trans_matr(stg_table,transition_rates_table_optim,''),...
	stg_sorting_cell,transition_rates_table_optim,x0),'x0');

[fcn_statsol_sum_sq_dev,~]=fcn_handles_fitting(y_data,x0,stg_table,stg_sorting_cell,nodes,predictor_names);

% default values for fitting hyperparameters:
% struct('CoolSched',@(T) (0.8*T), 'Generator',@(x) (x+(randperm(length(x))==length(x))*randn/100),...
% 'InitTemp',1,'MaxConsRej',1000, 'MaxSuccess',20,...
% 'MaxTries',300, 'StopTemp',1e-8, 'StopVal',-Inf, 'Verbosity',1);

% initial guess for parameters
init_par_vals=data_param_vals.*normrnd(1,1,size(predictor_names)); 
init_error=fcn_statsol_sum_sq_dev(init_par_vals);

% initial value of model nodes (with the initial parameter guess)
y_init=fcn_calc_init_stat_nodevals(x0,split_calc_inverse(fcn_build_trans_matr(stg_table,...
		fcn_trans_rates_table(nodes,'uniform',[],[],predictor_names,init_par_vals),''),...
		stg_sorting_cell,transition_rates_table_optim,x0),'');
% FIT
tic; [optim_par_vals,best_error,T_loss]=anneal(fcn_statsol_sum_sq_dev,init_par_vals,fitting_arguments); toc;

% RESULTS
% model variable values with fitted parameters
y_optim_param=fcn_calc_init_stat_nodevals(x0,...
  split_calc_inverse(fcn_build_trans_matr(stg_table,...
  fcn_trans_rates_table(nodes,'uniform',[],[],predictor_names,optim_par_vals),''),...
  stg_sorting_cell,transition_rates_table,x0),'');

% model variables: initial guess, true values (data), fitted values
data_init_optim=[y_init; y_data; y_optim_param]; 
min_val=min(min(data_init_optim(:,3:end))); max_val=max(max(data_init_optim(:,3:end)));
% parameters: initial guess, true values, fitted values
param_sets=[init_par_vals;data_param_vals;optim_par_vals];

%% PLOT

figure('name','param fitting (simul.ann.)'); 
% select nodes to plot (here we selected nodes that are not always 0 or 1)
sel_nodes=find(sum(data_init_optim)>0 & sum(data_init_optim)<3);
% PLOT fitting process
thres_ind=size(T_loss,1); % thres_ind=find(T_loss(:,2)<1e-2,1); 
plot_settings=[24 30];
figure('name','simul anneal')
fcn_plot_paramfitting(data_init_optim,T_loss,nodes,sel_nodes,[1 2],thres_ind,plot_settings)

%% Fitting by initial numerical gradient

[~,~,predictor_names]=fcn_get_trans_rates_tbl_inds(scan_params_sensit,scan_params_up_down_sensit,nodes); 
% create simulated data by randomly generating parameter values
data_param_vals=lognrnd(1,1,1,numel(predictor_names)); 
transition_rates_table_optim=fcn_trans_rates_table(nodes,'uniform',[],[],predictor_names,data_param_vals);

% generate 'data' from parameters
y_data=fcn_calc_init_stat_nodevals(x0,split_calc_inverse(...
	fcn_build_trans_matr(stg_table,transition_rates_table_optim,''),stg_sorting_cell,...
	transition_rates_table_optim,x0),'x0');

% anonymous functions for fitting
[~,fcn_statsol_values]=fcn_handles_fitting(y_data,x0,stg_table,stg_sorting_cell,nodes,predictor_names);

% initial values for parameters, model variables and error
init_par_vals=data_param_vals.*lognrnd(0,2,size(predictor_names));
init_vals=fcn_statsol_values(init_par_vals); init_error=sum((y_data-init_vals).^2); 

error_thresh=0.1; 	% what % of initial error to stop?
step_thresh=[]; 	% what step # to stop? you can leave this empty 
% init_error_table: changes to initial error when increasing or decreasing parameter values
init_error_table=[]; % if we have it from previous fitting than feed it to fcn
% incr_resol_init: initial % change from initial param values to calculate the numerical gradient
% incr_resol: change in param values during gradient descent
incr_resol_init=0.15; incr_resol=0.03;

% FIT
[init_error_table,optim_pars_conv,statsol_parscan,error_conv]=fcn_num_grad_descent(init_error_table,...
	{y_data,x0,stg_table,stg_sorting_cell,nodes,predictor_names},data_param_vals,...
	init_par_vals,incr_resol,incr_resol_init,error_thresh,[]);

%% PLOT

figure('name','numer grad_desc')
data_init_optim=[statsol_parscan([1 end],:); y_data];
fcn_plot_paramfitting(data_init_optim,error_conv,nodes,sel_nodes,[],[],plot_settings)

%% SAVE
fig_name=strcat('grad_descent',num2str(numel(predictor_names)),'fittingpars');
fcn_save_fig(fig_name,plot_save_folder,fig_file_type{3},'overwrite',resolution_dpi)
