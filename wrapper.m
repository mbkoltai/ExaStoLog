% This file contains the commands to run the functions calculating the
% stationary solution of stoch logical models, plot results and perform parametric analysis

path_to_toolbox='/bioinfo/users/mkoltai/research/models/KRAS_DNA_repair_model/matlab_ode/maboss_analytical/';
cd(path_to_toolbox)
% ADD FUNCTIONS to PATH
addpath('functions/') 
% using these publicly available MATLAB toolboxes:
% needed for PLOTS:
addpath('heatmaps') % https://mathworks.com/matlabcentral/fileexchange/24253-customizable-heat-maps
addpath('redblue'); % https://fr.mathworks.com/matlabcentral/fileexchange/25536-red-blue-colormap (for redblue colormaps)
% optional for PLOTS and saving of PLOTS
% for subplots with smaller gaps:
addpath('tight_subplot') % https://mathworks.com/matlabcentral/fileexchange/27991-tight_subplot-nh-nw-gap-marg_h-marg_w
% export figures as EPS or PDF as they appear on plots:
addpath('altmany-export_fig-acfd348') % Optional. https://mathworks.com/matlabcentral/fileexchange/23629-export_fig 
% optional for PARAMETER FITTING by simulated annealing:
addpath('anneal') % https://fr.mathworks.com/matlabcentral/fileexchange/10548-general-simulated-annealing-algorithm

set(0,'DefaultAxesTitleFontWeight','normal');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% model set-up

% models can be defined 
% A) by entering the list of nodes and their
% corresponding rules as a cell of strings, using MATLAB logical notation ('&', '|', '~', '(', ')'),
% for instance:
% nodes = {'cc','kras', 'dna_dam', 'chek1', 'mk2', 'atm_atr', 'hr','cdc25b', 'g2m_trans', 'cell_death'};
% rules={'cc',...
% 'kras',...
% '(dna_dam | kras) & ~hr',...
% 'atm_atr',...
% 'atm_atr & kras',...
% 'dna_dam',...
% '(atm_atr  | hr) & ~cell_death',...
% '(cc|kras) & (~chek1 & ~mk2) & ~cell_death',...
% 'g2m_trans | cdc25b',...
% 'cell_death | (dna_dam & g2m_trans)'}; 

% name of the model
model_name='kras15vars';

% OR B) the model can be read in from an existing BOOLNET file
[nodes,rules]=fcn_bnet_readin('krasmodel15vars.bnet'); % krasmodel10vars.bnet

% once we have the list of nodes and their logical rules, we can check if
% all variables referred to by rules are found in the list of nodes:
fcn_nodes_rules_cmp(nodes,rules)

% if yes, we generate a function file, which will create the model
truth_table_filename='fcn_truthtable.m'; fcn_write_logicrules(nodes,rules,truth_table_filename)

% from the model we generate the STG table, that is independent of values of transition rates
tic; [stg_table,~,~]=fcn_build_stg_table(truth_table_filename,nodes,'',''); toc
% [state_transitions_inds,K_sparse,A_sparse_fast]=fcn_build_stg_table(truth_table_filename,nodes,transition_rates_table,'num_matrix');

%% generate from already existing STG table (filling in values)

% to define transition rates, we can select given rates to have different values than 1, or from randomly chosen
% name of rates: 'u_nodename' or 'd_nodename'
chosen_rates={'u_cdc25b','d_dna_dam'}; chosen_rates_vals=[0.25, 0.15]; 

% then we generate the table of transition rates: first row is the 'up'
% rates, second row 'down' rates, in the order of 'nodes'
% ARGUMENTS
distr_type={'uniform','random'}; % uniform assigns a value of 1 to all params. other option: 'random'
meanval=[]; sd_val=[]; % if 'random' is chosen, the mean and standard dev of a normal distrib has to be defined 
% transition_rates_table=fcn_trans_rates_table(nodes,uniform_or_rand,meanval,sd_val,chosen_rates,chosen_rates_vals);
transition_rates_table=fcn_trans_rates_table(nodes,distr_type{1},meanval,sd_val,[],[]);
% meanval=1; sd_val=1; transition_rates_table=fcn_trans_rates_table(nodes,'random',meanval,sd_val,chosen_rates,chosen_rates_vals)

% build transition matrix A with parameter values
tic; [A_sparse,~]=fcn_build_trans_matr(stg_table,transition_rates_table,''); toc
% if we want the kinetic matrix too, this is the 2nd output of the function
% tic; [A_sparse,K_sparse]=fcn_build_trans_matr(stg_table,transition_rates_table,'kinetic'); toc

% density of transition matrix A
% nnz(A_sparse)/numel(A_sparse)

%% calculate steady states with kernel/nullspace calculation functions

% defining an initial condition
n_nodes=numel(nodes); truth_table_inputs=rem(floor([0:((2^n_nodes)-1)].'*pow2(0:-1:-n_nodes+1)),2);
% define initial values
x0=zeros(1,2^n_nodes)'; 
% defining a dominant initial state (eg. dom_prob=0.8, ie. 80% probability)
dom_prob=0.8; initial_state=[1 1 1 zeros(1,numel(nodes)-3)];
x0(ismember(truth_table_inputs,initial_state,'rows'))=dom_prob; 
% take those states where first 3 variables have a value of 1, and we want them to have a nonzero probability
sel_states=all(truth_table_inputs(:,1:3)');
% create a vector of random probabilities for these states, with a sum of (1-dom_prob)
rand_vect=abs(rand(sum(sel_states)-1,1)); rand_vect=rand_vect/( sum(rand_vect)/(1-dom_prob) );
x0(~ismember(truth_table_inputs,initial_state,'rows') & sel_states') = rand_vect;

% completely random initial condition
% x0=rand(1,size(truth_table_inputs,1))'; x0=x0/sum(x0);

% CALCULATE STATIONARY STATE
% ARGUMENTS:
% transition matrix: A
% table of transition rates: transition_rates_table
% initial conditions: x0
tic; [stat_sol,term_verts_cell,cell_subgraphs]=split_calc_inverse(A_sparse,transition_rates_table,x0); toc
% OUTPUTS
% stat_sol: stationary solution for all the states
% term_verts_cell: index of nonzero states. If the STG is disconnected, then the nonzero states corresp to these disconn subgraphs are in separate cells
% cell_subgraphs: indices of states belonging to disconnected subgraphs (if any)

% nonzero states can be quickly queried by:
stat_sol(stat_sol>0) % probability values of nonzero states
truth_table_inputs(stat_sol>0,:) % logical states that are nonzero

% sum the probabilities of nonzero states by nodes, both for the initial condition and the stationary solution
% ARGUMENTS
% initial conditions: x0
% stat_sol: stationary solution for all the states
% nodes: list of nodes
[stationary_node_vals,init_node_vals]=fcn_calc_init_stat_nodevals(x0,stat_sol);

% can check results by doing direct matrix exponentiation, for systems up to ~10 nodes, but only if there are no cycles.
% x_sol=((x0')*A_sparse^1e5)';
% stationary_node_vals=x_sol'*truth_table_inputs;
% checked with MaBoSS simuls, results are identical (up to 1% dev.) as they have to be

%% PLOTTING RESULTS

% set(0,'DefaultFigureWindowStyle','docked')

% PLOT A/K and stat solutions
% ARGUMENTS
% K_sparse or A_sparse: kinetic or transition matrix
% min_max_col: minimum and max color for heatmap
% fontsize: [font size on the heatmap, title font size for stationary solutions]
% barwidth_states_val: width of the bars for bar plot of stationary solutions of states
% sel_nodes: nodes to show. If left empty, all nodes are shown
sel_nodes=[]; min_max_col=[-1 1];barwidth_states_val=2;fontsize=[10 20]; % fontsize_hm,fontsize_stat_sol
% CAREFUL! if more than 12 nodes, generating the figure for A_sparse can be time-consuming
fcn_plot_A_K_stat_sol(A_sparse, nodes, sel_nodes, stat_sol, x0, min_max_col,fontsize,barwidth_states_val,[])

% PLOT stationary solutions (without A/K matrix)
% nonzero_flag: if non-empty, only the nonzero states are shown
sel_nodes=3:numel(nodes); nonzero_flag='nonzero_only';
barwidth_states_val=0.8; % for 3 nonzero states ~0.8 is a good value
% if init_node_vals='', then only shows stationary solution
% nonzero_flag: if this is non-empty, then we only plot the nonzero states, this is useful for visibility if there are many states
fcn_plot_A_K_stat_sol([], nodes, sel_nodes, stat_sol, x0, min_max_col,fontsize,barwidth_states_val,nonzero_flag)

% SAVE
% export_fig kras_stationary_sols.eps -transparent

% PLOT binary heatmap of nonzero stationary states by NODES
% ARGUMENT
% term_verts_cell: which subgraph to plot if there are disconnected ~
% num_size_plot: font size of 0/1s on the heatmap
% hor_gap: horizontal gap between terminal SCCs, bottom_marg: bottom margin, left_marg: left margin
numsize_plot=16; fontsize=12; hor_gap=0.02; bottom_marg=0.14; left_marg=0.06; 
param_settings=[numsize_plot fontsize hor_gap bottom_marg left_marg];
% idnex of nonempty subgraph, check by <term_verts_cell>
nonempty_subraph=4;
% want to use tight subplot?
tight_subplot_flag='yes';
% nodes to show
sel_nodes=3:numel(nodes);
% probability threshold for states to show (if left empty, all states shown)
prob_thresh=0.05;  % []; % 0.05;
% PLOT
statsol_binary_heatmap=fcn_plot_statsol_bin_hmap(stat_sol,prob_thresh,term_verts_cell{nonempty_subraph},nodes,sel_nodes,param_settings,tight_subplot_flag);
% export_fig sink_states_kras_model_heatmap.eps -transparent

%% plot of STG (state transition graph)

% for models larger than 12 nodes generating these plots can be very time consuming

% PLOT the STG of the model. 
% show the full STG and a selected (counter) subgraph
subgraph_index=find(cellfun(@(x) numel(x), term_verts_cell)>0); % select non-empty subgraph
titles = {'Full state transition graph',strcat('subgraph #',num2str(subgraph_index))}; 
% cropping the subplots (optional)
xlim_vals=[0 21;-5 5]; ylim_vals=[0 23;-5 5]; 
% parameters for plot
default_settings=[20 1 7 5 8]; % fontsize, linewidth_val, arrowsize, default_markersize, highlight_markersize
% color of source states of STG
source_color='blue'; 
% figure(); 
plot_STG(A_sparse,subgraph_index,default_settings,xlim_vals,ylim_vals,titles,source_color)

% PLOT a single STG (that can be a selected subgraph of entire STG)
% cropping (optional)
xlim_vals=[-4 5]; ylim_vals = [-5 5]; 
titles ={strcat('subgraph #',num2str(subgraph_index))};
A_sub=A_sparse(cell_subgraphs{subgraph_index},cell_subgraphs{subgraph_index});
default_settings=[20 1 7 5 8]; % fontsize,linewidth_val, arrowsize, default_markersize, highlight_markersize
plot_STG(A_sub,'',default_settings,xlim_vals,ylim_vals,titles,source_color)
% plot_STG(A_sub,'',default_settings,[],[],titles,source_color)

%% STGs on subplots, with given parameter highlighted on each

% if STG table doesn't exist generate with
% [stg_table,~,~]=fcn_build_stg_table(truth_table_filename,nodes,'','');

selected_pars=[1 3 4 5 6 11 9 8]; % parameters to highlight, either 'all', or numeric array [1 2 3]
plot_pars=[20 0.1 7 3 6]; % plot_pars=[fontsize,linewidth_val, arrowsize, default_markersize, highlight_markersize]
% parameters for highlighted transitions: color and width of corresp edges
highlight_settings={'yellow',3}; 
% if using tight_subplots toolbox:
% [ha,~] = tight_subplot(4,4,[0.06 0.02],[0.05 0.05],[0.05 0.05]);
% tight_subplot(Nh, Nw, gap, marg_h, marg_w): 
% gap: gaps between the axes in normalized units
% marg_h: margins in height in normalized units
% marg_w: margins in width in normalized units
tight_subpl_flag='yes'; tight_subplot_pars=[0.06 0.02; 0.05 0.05; 0.05 0.05]; 
% cropping plot (optional, for better visibility)
limits=[-4 5;-5 5]; 

figure()
plot_STG_sel_param(A_sparse,counter,nodes,cell_subgraphs,selected_pars,stg_table,plot_pars,highlight_settings,limits,tight_subpl_flag,tight_subplot_pars)
% show all params that have an effect
% plot_STG_sel_param(A_sparse,counter,nodes,cell_subgraphs,'all',stg_table,plot_pars,highlight_settings,'',tight_subpl_flag,tight_subplot_pars)
% SAVE
% export_fig kras_model_krasMUT_STG_colorededges_thick.eps -transparent -nocrop

% only one graph: in this case this is a subgraph of the entire STG, but could be the entire STG too
counter=4; B=A_sparse(cell_subgraphs{counter},cell_subgraphs{counter});
B_state_transitions_inds=fcn_stg_table_subgraph(stg_table,cell_subgraphs,counter);
% 2nd, 4th arguments left empty
% figure(); 
plot_STG_sel_param(B,'',nodes,'','all',B_state_transitions_inds,default_settings,highlight_settings,limits,'yes',tight_subplot_pars)

%% Parameter sensitivity analysis: one-dimensional parameter scans

% select the nodes whose parameters we want to scan in:
% all nodes that have actual transitions
par_inds_table=unique(stg_table(:,[3 4]),'rows');
scan_params=unique(par_inds_table(:,1))'; 
% or by selecting given nodes by their names: scan_params=find(ismember(nodes,{'cdc25b','atm_atr'}));
scan_params_up_down=arrayfun(@(x) par_inds_table(par_inds_table(:,1)==x,2)', scan_params,'un',0); 
% num2cell(repelem([1 2],numel(scan_params),1),2)'; % both up and down rates
% num2cell(ones(1,numel(scan_params))); % only up 
% num2cell(2*ones(1,numel(scan_params))); % only down
% {[1 2], 1, [1 2]}; % manually selected

% min and max of range of values; resolution of the scan; linear or logarithmic sampling
parscan_min_max = [1e-2 1e2]; n_steps=10; sampling_types={'log','linear'}; 

% FUNCTION for generating matrix of ordered values for the parameters to scan in
% [scan_par_table,scan_par_inds,~]= fcn_get_trans_rates_tbl_inds(scan_params,scan_params_up_down,transition_rates_table);
parscan_matrix=fcn_onedim_parscan_generate_matrix(scan_params,scan_params_up_down,nodes,sampling_types{1},parscan_min_max,n_steps);
% set values in one/more columns manually (order is {3,1}->{3,2}->{4,1}->{4,2} etc)
% parscan_matrix(:,2)=logspace(-1,1,size(parscan_matrix,1));
% entire matrix can be made a random matrix

% FUNCTION for 1-DIMENSIONAL PARSCAN
% ARGUMENTS
% stg table: generated above, state transition graph
% transition_rates_table:default values of transition rates
% initial conditions: x0
% stg_table: generated by [stg_table,~,~]=fcn_build_stg_table(truth_table_filename,nodes,transition_rates_table,'');
% [~,scan_par_inds,~]=fcn_get_trans_rates_tbl_inds(scan_params,scan_params_up_down,nodes);
[stationary_state_vals_onedimscan,stationary_node_vals_onedimscan,stationary_state_inds_scan]=fcn_onedim_parscan_calc(stg_table,transition_rates_table,...
                                                            x0,nodes,parscan_matrix,scan_params,scan_params_up_down);

%% PLOT RESULTS of 1-by-1 parameter scan on heatmap/lineplot

%%% SECOND PLOT TYPE: show the stationary value or response coefficient of 1 variable or state on 1 subplot, as a fcn of all relevant parameters
nonzero_states_inds=find(stat_sol>0);
sensit_cutoff=0.5; % minimal value for response coefficient (local sensitivity)
% nonzero states of the model
nonzero_states=unique(cell2mat(stationary_state_inds_scan(:)'))';
% select parameters of plot
height_width_gap=[0.1 0.03]; bott_top_marg =[0.05 0.1]; left_right_marg=[0.05 0.05];
% [fontsize_axes,fontsize_title,params_tight_subplots(leave empty if not installed),model_name]
plot_param_settings={12,14,{height_width_gap bott_top_marg left_right_marg},model_name}; % plot_param_settings={12,14,[],model_name}; 
% select type of plot
plot_types={{'lineplot','heatmap'} {'nodes','states'} {'values','sensitivity'}};
all_opts_perm=[[1 1 1]; unique([perms([1 1 2]); perms([2 2 1])],'rows'); [2 2 2]];
for k=1:size(all_opts_perm,1)
plot_type_options=all_opts_perm(k,:);
figure
[resp_coeff,scan_pars_sensit,scan_params_sensit_up_down,fig_filename]=fcn_onedim_parscan_plot_parsensit(plot_types,plot_type_options,...
                                                           stationary_node_vals_onedimscan,stationary_state_vals_onedimscan,...
                                                           nonzero_states_inds,parscan_matrix,nodes,scan_params,scan_params_up_down,...
                                                           sensit_cutoff,plot_param_settings);

% <resp_coeffs> dimensions: (parameters, values,nodes), so resp_coeffs(:,:,7) are the resp. coeff values across the param ranges for the 7th node

save_folder='sample_plots/'; fig_file_type='.png'; % '.eps'
export_fig(strcat(save_folder,fig_filename,fig_file_type),'-transparent','-nocrop')
end

%% multidimensional parameter scan: LATIN HYPERCUBE SAMPLING (random multidimensional sampling within given parameter ranges)

% WHICH PARAMETERS to scan it simultaneously?
% all relevant parameters
% scan_params=unique(stg_table(:,3))'; scan_params_up_down=num2cell(repelem([1 2],numel(scan_params),1),2)'; 
% scan_params: scanned parameters, id by node. 
scan_params=scan_pars_sensit; % scan_params=unique(stg_table(:,3))';
% scan_params_up_down: id by up (1) or down (2) rate
scan_params_up_down=scan_params_sensit_up_down; % num2cell(ones(1,numel(scan_params))) {[1 2], 1, [1 2]}; 
% scan_params_sensit_up_down=arrayfun(@(x) scan_params_sensit_up_down{x}', 1:numel(scan_params_sensit_up_down),'un',0)

% PERFORM Latin Hypercube sampling (LHS) SAMPLING
sampling_types={'lognorm','linear','logunif'};
sampling_type=sampling_types{3};
% <lhs_scan_dim>: number of param sets
lhs_scan_dim=1e3;
% par_min_mean: minimum or in case of lognormal the mean of distribution. Can be a scalar or a vector, if we want different values for different parameters
% max_stdev: maximum or in case of lognormal the mean of distribution. Can be a scalar or a vector, if we want different values for different parameters
% for 'lognorm' and 'logunif' provide the log10 value of desired mean/min and stdev/max!!
par_min_mean=-2; % repmat(1.5,1,numel(cell2mat(scan_params_up_down(:)))); par_min_mean(4)=3; 
max_stdev=2; % repmat(0.5,1,numel(cell2mat(scan_params_up_down(:))));
tic;
[all_par_vals_lhs,stat_sol_lhs_parscan,...
    stat_sol_states_lhs_parscan,stat_sol_states_lhs_parscan_cell]=fcn_multidim_parscan_latinhypcube(par_min_mean,max_stdev,sampling_type,...
                                                    lhs_scan_dim,scan_params,scan_params_up_down,transition_rates_table,stg_table,x0,nodes);
toc;
% <all_par_vals_lhs>: table of parameter sets
% <stat_sol_lhs_parscan>: stationary values of nodes
% <stat_sol_states_lhs_parscan>: stationary values of states
% <stat_sol_states_lhs_parscan_cell>: indices of non-zero states, if they are always the same the variable is a 1-row array, 
% otherwise a cell with the non-empty elements stat_sol_states_lhs_parscan_cell{k,1} contain values of nonzero
% states, elements stat_sol_states_lhs_parscan_cell{k,2} their indices 

%% SCATTERPLOTS of STATE or NODE values as a function of the selected parameters, with the trendline shown (average value per parameter bin)

var_ind=7; % which STATE or NODE to plot
% <all_par_vals_lhs>: parameter sets
param_settings = [50 4 16 stat_sol_states_lhs_parscan_cell]; % [number_bins_for_mean,trendline_width,axes_fontsize,index nonzero states]
% STATES or NODES? - <scan_values>: values to be plotted
scan_values=stat_sol_lhs_parscan; % stat_sol_states_lhs_parscan
% PLOT
sampling_type=sampling_types{3}; % sampling_types={'lognorm','linear','logunif'};
fcn_multidim_parscan_scatterplot(var_ind,all_par_vals_lhs,scan_values,scan_params,scan_params_up_down,nodes,sampling_type,param_settings) 

save_folder='sample_plots/'; fig_file_type='.png'; % '.eps'
fig_name=strcat(save_folder,'multidim_parscan_scatterplot_',nodes{var_ind},'_',model_name,fig_file_type);
export_fig(fig_name,'-transparent','-nocrop')

%% calculating & plotting (heatmap) correlations between variables OR between params and variables by linear/logist regression

% ARGUMENTS
% stat_sol_lhs_parscan: values from parameter sampling
% nodes: name of ALL nodes
% sel_nodes: name of selected nodes (pls provide in ascending order)
% fontsize: ~ for labels and titles (displaying correlation)
% HEATMAPS of correlations between selected variables
sel_nodes=3:15; plot_settings=[15 16]; % [fontsize on plot, fontsize on axes/labels]
plot_type_flag={'var_var','heatmap'}; % this is plotting the heatmap of correlations between variables
[varvar_corr_matr,p_matrix_vars]=fcn_multidim_parscan_parvarcorrs(plot_type_flag,all_par_vals_lhs,stat_sol_lhs_parscan,...
                                            nodes,sel_nodes,scan_params,scan_params_up_down,[],plot_settings);

save_folder='sample_plots/'; fig_file_type={'.png','.eps'};
fig_name=strcat(save_folder,model_name,'_',strjoin(plot_type_flag,'_'),fig_file_type{2});
export_fig(fig_name,'-transparent','-nocrop')
                   
%% scatterplots of selected variables [i,j]: var_i VS var_j

sel_nodes=10:15; plot_settings=[10 12]; % [fontsize_axes, fontsize_titles]
plot_type_flag={'var_var','scatter'}; % this is plotting the scatterplots of variables with correlation values
fcn_multidim_parscan_parvarcorrs(plot_type_flag,all_par_vals_lhs,stat_sol_lhs_parscan,...
                                    nodes,sel_nodes,scan_params,scan_params_up_down,[],plot_settings);

save_folder='sample_plots/'; fig_file_type={'.png','.eps'};
fig_name=strcat(save_folder,model_name,'_',strjoin(plot_type_flag,'_'),fig_file_type{2});
export_fig(fig_name,'-transparent','-nocrop')
                                
%% linear or lin-log regression of VARIABLES as fcn of PARAMETERS: VARIABLE=f(PARAMETER), the function plots R squared

plot_type_flag={'par_var','heatmap','r_sq'}; % {'par_var','heatmap'/'lineplot','r_sq'/'slope'}
sel_nodes=setdiff(3:numel(nodes),8); 
% plot_settings=[fontsize,maximum value for heatmap colors], if plot_settings(2)=NaN, then max color automatically selected
plot_settings=[14 NaN]; 
% if regression type is 'linlog', then the fit is y = a + b*log10(x)
regr_type={'log','linear'}; % linlog recommended if parameter values log-uniformly distributed in sampling
[r_squared,slope_intercept]=fcn_multidim_parscan_parvarcorrs(plot_type_flag,all_par_vals_lhs,stat_sol_lhs_parscan,...
                                        nodes,sel_nodes,scan_params,scan_params_up_down,regr_type{1},plot_settings);

save_folder='sample_plots/'; fig_file_type={'.png','.eps'};
fig_name=strcat(save_folder,model_name,'_',strjoin(plot_type_flag,'_'),fig_file_type{1});
export_fig(fig_name,'-transparent','-nocrop')
                                    
%% Quantify importance of parameters from LHS by a regression tree

% predictor importance values: look into arguments of <fitrtree> in MATLAB documentation to modulate regression tree
% for STATES or NODES?
scan_values=stat_sol_lhs_parscan; % stat_sol_states_lhs_parscan
sel_nodes=3:numel(nodes); % STATES or NODES to be analyzed
% names of selected transition rates and their predictor importance values
% [~,~,predictor_names] = fcn_get_trans_rates_tbl_inds(scan_params,scan_params_up_down,nodes);
% predictorImportance_vals=cell2mat(arrayfun(@(x) predictorImportance(...
%     fitrtree(all_par_vals_lhs,scan_values(:,x),'PredictorNames',predictor_names)), sel_nodes,'un',0)');

% CALCULATE and PLOT predictor importance
plot_type_flags={'line','bar'};
[predictor_names,predictorImportance_vals]=fcn_multidim_parscan_predictorimport(scan_params,scan_params_up_down,...
                                                all_par_vals_lhs,scan_values,nodes,sel_nodes,plot_type_flags{2});

save_folder='sample_plots/'; fig_file_type={'.png','.eps'};
fig_name=strcat(save_folder,model_name,'_','regression_tree_pred_import',fig_file_type{2});
export_fig(fig_name,'-transparent','-nocrop')
                                            
%% Sobol total sensitivity metric                                            

% On Sobol total sensitivity index see: https://en.wikipedia.org/wiki/Variance-based_sensitivity_analysis
% This metric indicates how much of the total variance in a variable is due to variation in a given parameter
% We calculate here the usual numerical approximation of analytical equivalent from Monte Carlo sampling.
% From the LHS sampling above we take the matrices of parameter sets and variable values:
% [parameter sets, variable values]: [all_par_vals_lhs,stat_sol_lhs_parscan]

% Sobol total sensitivity: calculated for one variable at a time
sel_vars=setdiff(1:numel(nodes),find(sum(cell2mat(arrayfun(@(x) strcmp(nodes,x), {'cc','KRAS','CDC25B'},'un',0)')))); % selected nodes to display
% sel_vars=[]; % if left empty, all nodes/states are analyzed
sample_size=500; % if left empty, the sample size is half of the original param scan <all_par_vals_lhs>
plot_settings=[14 14 22 NaN NaN 10]; % [fontsize_plot,fontsize_axes,fontsize_title, min_color(optional), max_color(opt), progress_calcul_every_x_% (opt)];
var_types={'nodes','states'}; % analysis for states or nodes
% to calculate Sobol total sensitivity we need <sample_size*numel(scan_params_up_down)> evaluations of the model
tic; sobol_sensit_index=fcn_multidim_parscan_sobol_sensit_index([],var_types{1},all_par_vals_lhs,stat_sol_lhs_parscan,stat_sol_states_lhs_parscan,sample_size,...
                                scan_params,scan_params_up_down,stg_table,x0,nodes,sel_vars,plot_settings); toc;

% if we have already calculated <sobol_sensit_index> and only want to plot results, provide it as the FIRST argument <sobol_sensit_index>
fcn_multidim_parscan_sobol_sensit_index(sobol_sensit_index,var_types{1},[],[],[],[],...
                                scan_params,scan_params_up_down,[],[],nodes,sel_nodes,plot_settings);

save_folder='sample_plots/'; fig_file_type={'.png','.eps'};
fig_name=strcat(save_folder,model_name,'_','sobol_sensitivity_index',fig_file_type{1});
export_fig(fig_name,'-transparent','-nocrop')
                            
%% PARAMETER FITTING

% recall relevant functions:
% transition_rates_table=fcn_trans_rates_table(nodes,'uniform',[],[],chosen_rates,chosen_rates_vals); % transition_rates_table=ones(size(transition_rates_table));
% tic; [A_sparse,~]=fcn_build_trans_matr(stg_table,transition_rates_table,''); toc; 
% tic; [stat_sol,term_verts_cell,cell_subgraphs]=split_calc_inverse(A_sparse,transition_rates_table,x0); toc
% [stationary_node_vals,init_node_vals]=fcn_calc_init_stat_nodevals(x0,stat_sol);

% define parameters to vary (predictor_names)
[~,~,predictor_names]=fcn_get_trans_rates_tbl_inds(scan_params,scan_params_up_down,nodes);
% define data vector (generate some data OR load from elsewhere)
sel_param_vals=lognrnd(1,1,1,numel(predictor_names)); % abs(normrnd(1,0.5,1,numel(predictor_names)));
transition_rates_table=fcn_trans_rates_table(nodes,'uniform',[],[],predictor_names,sel_param_vals);
y_data=fcn_calc_init_stat_nodevals(x0,split_calc_inverse(fcn_build_trans_matr(stg_table,transition_rates_table,''),transition_rates_table,x0));

% create functions that calculate sum of squared deviations & values of
% variables (composed of different fcns) - RERUN THIS if you want to fit to new data or new non-fitted transition rates!!
[fcn_statsol_sum_sq_dev,fcn_statsol_values]=fcn_simul_anneal(y_data,x0,stg_table,nodes,predictor_names);
% evaluate/test: abs(fcn_statsol_sum_sq_dev(chosen_rates_vals) - fcn_statsol_sum_sq_dev(zeros(size(chosen_rates_vals))) )

% FITTING by simulated annealing (look at arguments in anneal/anneal.m)
% initial guess for parameters
init_vals=rand(size(predictor_names)); init_error=fcn_statsol_sum_sq_dev(init_vals);
% simulated annealing with existing algorithm anneal/anneal.m, with modifications in script: defined counter=0 before while loop, 
% and inserted <T_loss(counter,:)=[T oldenergy];> at line 175, defined <T_loss> as 3rd output
tic; [optim_par_vals,best_error,T_loss]=anneal(fcn_statsol_sum_sq_dev,init_vals,struct('Verbosity',2)); toc % 'StopTemp',1e-8
% PLOT results
thres=1e-4; semilogy(1:find(T_loss(:,2)<thres,1), T_loss(1:find(T_loss(:,2)<thres,1),:),'LineWidth',4); legend({'temperature', 'SSE'},'FontSize',22);
xlabel('number of iterations','FontSize',16); set(gca,'FontSize',16); title('Parameter fitting by simulated annealing','FontSize',22)

% normalized by values: mean absolute fractional error
mean( abs(y_data - fcn_statsol_values(optim_par_vals))./y_data )
% distance of found params from true values (fractional difference)
abs(optim_par_vals - sel_param_vals)./sel_param_vals