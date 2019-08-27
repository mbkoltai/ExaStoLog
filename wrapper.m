% This file contains the commands to run the functions calculating the
% stationary solution of stoch logical models, plot results and perform parametric analysis

% go to the folder of the file
editor_service = com.mathworks.mlservices.MLEditorServices; editor_app = editor_service.getEditorApplication;
active_editor = editor_app.getActiveEditor; storage_location = active_editor.getStorageLocation;
file = char(storage_location.getFile); path_to_toolbox = fileparts(file); cd(path_to_toolbox);

% ADD FUNCTIONS to PATH
add_functions

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% model set-up

% models can be defined 
% A) by entering the list of nodes and their
% corresponding rules as a cell of strings, using MATLAB logical notation ('&', '|', '~', '(', ')'),
% for instance a toy model with cyclic attractor:
% nodes={'A','B','C'}; % rules={'~A','A','A&~C'}
% nodes={'A','B','C','D'}; rules={'~B','~A&C','B|C&~D','C'}

% LIST of MODELS
model_name_list = {'mammalian_cc', ...
'krasmodel15vars', ...
'breast_cancer_zanudo2017'....
'dnarepair_rodriguez_15nodes',...
'EMT_cohen_ModNet'}; % 
% name of the model
model_index=5;
model_name=model_name_list{model_index};

% where to save figures
plot_save_folder=strcat('doc/sample_plots/',model_name,'/');

% model read in from an existing BOOLNET file
[nodes,rules]=fcn_bnet_readin(strcat('model_files/',model_name,'.bnet')); 

% once we have list of nodes and their logical rules, check if all variables referred to by rules found in list of nodes:
fcn_nodes_rules_cmp(nodes,rules)

% if yes, we generate a function file, which will create the model
truth_table_filename='fcn_truthtable.m'; fcn_write_logicrules(nodes,rules,truth_table_filename)

% from the model we generate the STG table, that is independent of values of transition rates
% this takes 20-30 seconds for 20 node model
tic; stg_table=fcn_build_stg_table(truth_table_filename,nodes); toc

% density of transition matrix
size(stg_table,1)/(2^(2*numel(nodes)))
% visualize transition matrix
% spy(A_sparse); xlabel('model states'); ylabel('model states'); set(gca,'FontSize',24)
% save: export_fig(strcat(save_folder,model_name,'_A_sparse.pdf'),'-transparent','-nocrop','-r350')

%% generate transition matrix from existing STG table

% to define transition rates, we can select given rates to have different values than 1, or from randomly chosen
% name of rates: 'u_nodename' or 'd_nodename'
% for KRAS model, state of KRAS defines whether its mutant or not: chosen_rates={'d_KRAS','d_cc'}; chosen_rates_vals=[0 0]; 
chosen_rates=[]; chosen_rates_vals=[];
% OR leave them empty: chosen_rates=[]; chosen_rates_vals=[];

% then we generate the table of transition rates: first row is the 'up'rates, second row 'down' rates, in the order of 'nodes'
% ARGUMENTS
distr_type={'uniform','random'}; % <uniform> assigns a value of 1 to all params. other option: <random>
meanval=[]; sd_val=[]; % if 'random' is chosen, the mean and standard dev of a normal distrib has to be defined 
transition_rates_table=fcn_trans_rates_table(nodes,distr_type{1},meanval,sd_val,chosen_rates,chosen_rates_vals);
% meanval=1; sd_val=1; transition_rates_table=fcn_trans_rates_table(nodes,'random',meanval,sd_val,chosen_rates,chosen_rates_vals)

% build transition matrix A with parameter values
tic; [A_sparse,~]=fcn_build_trans_matr(stg_table,transition_rates_table,''); toc
% if we want the kinetic matrix too, this is the 2nd output of the function
% tic; [A_sparse,K_sparse]=fcn_build_trans_matr(stg_table,transition_rates_table,'kinetic'); toc
%

% VISUALIZE transition matrix
% spy(A_sparse); xlabel('model states'); ylabel('model states'); set(gca,'FontSize',24)

%% define initial condition

n_nodes=numel(nodes); 
% truth_table_inputs=rem(floor([0:((2^n_nodes)-1)].'*pow2(0:-1:-n_nodes+1)),2);

% define some nodes with a fixed value and a probability <dom_prob>: 
% states satisfying this condition will have a total initial probability of <dom_prob>, the other states 1-dom_prob

% CELL CYCLE MODEL
% initial state specifying all variables:
% initial_fixed_nodes={'CycD','Rb_b1','Rb_b2','p27_b1','p27_b2','Cdh1','Skp2','E2F','CycE','CycA','CycB','Cdc20','UbcH10'}; 
% initial_fixed_nodes_vals=[ones(1,7) zeros(1,6)];
% 
% initial state specifying some variables:: CycE=0 & CycA=0 & CycB=0 & Cdh1=1 & % Rb=1 & p27=1
% initial_fixed_nodes={'CycE','CycA','CycB','Cdh1','Rb_b1','Rb_b2','p27_b1','p27_b2'}; initial_fixed_nodes_vals=[0 0 0 1 1 1 1 1];
%
%
% KRAS model, WT: {'cc','KRAS'}, [1 0]. Mutant: {'cc','KRAS'}, [1 0]
% KRAS 15 nodes mutant
% initial_fixed_nodes={'cc','KRAS','cell_death'}; initial_fixed_nodes_vals=[1 1 0];

initial_fixed_nodes_list={ {'CycE','CycA','CycB','Cdh1','Rb_b1','Rb_b2','p27_b1','p27_b2'}, ... % mammalian_cc
                             {'cc','KRAS','DSB','cell_death'}, ...                              % krasmodel15vars
                              {'Alpelisib', 'Everolimus','PIM','Proliferation','Apoptosis'},...  % breast_cancer_zanudo2017 % ,'PIM'
                          {''},... % Rodriguez
                          {'ECMicroenv','DNAdamage','Metastasis','Migration','Invasion','EMT','Apoptosis'}}; % 
                          % {'Alpelisib', 'Everolimus', 'PI3K', 'PIM', 'PDK1', 'Proliferation', 'Apoptosis'}

initial_fixed_nodes_vals_list = {[0 0 0 1 1 1 1 1], ... % mammalian_cc
    [1 1 1 0], ... % krasmodel15vars: [1 1] is cell cycle ON, KRAS mutation ON
    [0 1 0 zeros(1,2)],...  % breast_cancer_zanudo2017
    [],...  % Rodriguez
    [1 1 zeros(1,5)]}; % EMT-Cohen model
initial_fixed_nodes=initial_fixed_nodes_list{model_index}; initial_fixed_nodes_vals=initial_fixed_nodes_vals_list{model_index};

% what is the probability of this state, (eg. dom_prob=0.8, ie. 80% probability)
dom_prob=1;
% if <random> the probability is randomly distributed among states, if <uniform> uniformly
distrib_types={'random','uniform'};
% if plot_flag non-empty, we get a bar plot of initial values
plot_flag='';
% function assigns a probability of <dom_prob> to the states with the fixed nodes having the defined values
tic; x0=fcn_define_initial_states(initial_fixed_nodes,initial_fixed_nodes_vals,dom_prob,nodes,distrib_types{1},plot_flag); toc

% completely random initial condition: 
% x0=zeros(2^n_nodes,1); x0=rand(size(truth_table_inputs,1),1); x0=x0/sum(x0);
% completely uniform initial condition
% x0=ones(2^numel(nodes),1)/(2^numel(nodes));

%% CALCULATE STATIONARY STATE

% ARGUMENTS:
% transition matrix: A
% table of transition rates: transition_rates_table
% initial conditions: x0
% get the subnetworks of the STG and topologically sort them
tic; stg_sorting_cell=fcn_scc_subgraphs(A_sparse,x0); toc
% <stg_sorting_cell> contains
% {subnetws: which subgraph each state belongs to,
% scc_submat_cell: states in the subgraphs,
% nonempty_subgraphs: which subgraphs are populated by the initial condition,
% sorted_vertices_cell: states (vertices) topologically sorted in each nonempty subgraph,
% cyclic_sorted_subgraphs_cell: sorted states within cycles}

% calculated stationary solution
tic; [stat_sol,term_verts_cell,cell_subgraphs]=split_calc_inverse(A_sparse,stg_sorting_cell,transition_rates_table,x0); toc
% OUTPUTS
% stat_sol: stationary solution for all the states
% term_verts_cell: index of nonzero states. If the STG is disconnected the nonzero states corresp to these disconn subgraphs are in separate cells
% cell_subgraphs: indices of states belonging to disconnected subgraphs (if any)

% nonzero states displayed:
stat_sol(stat_sol>0)' % probability values of nonzero states
truth_table_inputs(stat_sol>0,:) % list of logical states that are nonzero

% sum the probabilities of nonzero states by nodes, both for the initial condition and the stationary solution
% ARGUMENTS
% initial conditions: x0
% stat_sol: stationary solution for all the states
% nodes: list of nodes
[stationary_node_vals,init_node_vals]=fcn_calc_init_stat_nodevals(x0,stat_sol,'x0');

% Checked with MaBoSS simuls, results are identical (up to 1% dev.).
% Comparing with simulation of mammalian cell cycle model with 12 nodes: 
% look in folder 'doc/sample_plots/maboss'

%% PLOTTING RESULTS

% this function plots 2 or 3 subplots:
% 1) transition matrix (optional)
% 2) stationary probability of the model's states
% 3) stationary probability of the model's variables (having the value of 1)

% PLOT A/K and stat solutions
% ARGUMENTS
% matrix_input: [], K_sparse or A_sparse (kinetic or transition matrix)
% min_max_col: minimum and max color for heatmap
% fontsize: [font size on the heatmap, title font size for stationary solutions]
% barwidth_states_val: width of the bars for bar plot of stationary solutions of states
% sel_nodes: nodes to show. If left empty, all nodes are shown
% nonzero_flag: minimal value for probability to display - if this is non-empty, only plot nonzero states, useful for visibility if there are many states
sel_nodes=[];  % 3:numel(nodes)
min_max_col=[0 1]; barwidth_states_val=0.8;fontsize=[10 20]; % fontsize_hm,fontsize_stat_sol
plot_settings = [fontsize barwidth_states_val min_max_col]; prob_thresh=0.03;
% WARNING!!! if more than 12 nodes, generating the figure for A/K can be time-consuming
matrix_input=A_sparse; % leave this variable empty to have only 2 subplots, without transition matrix
figure('name','A_K_stat_sol')
fcn_plot_A_K_stat_sol(matrix_input,nodes,sel_nodes,stat_sol,x0,plot_settings,prob_thresh)

% SAVE
% enter any string for the last argument to overwrite existing plot!!
if exist(plot_save_folder,'dir')==0; mkdir(plot_save_folder); end
fig_file_type={'.png','.eps','.pdf','.jpg','.tif'}; 
if ~isempty(matrix_input);  matrix_input_str='_with_matrix'; else; matrix_input_str=''; end
% if <overwrite_flag> non-empty then existing file with same name is overwritten. 
overwrite_flag='yes'; 
% for <resolution> you can enter dpi value manually
resolution_dpi='-r350'; % magnification=0.8;  strcat('-r',num2str(magnification*get(0, 'ScreenPixelsPerInch')));
fcn_save_fig(strcat('single_solution_states_nodes_stat_sol',matrix_input_str),plot_save_folder,fig_file_type{3},overwrite_flag,resolution_dpi);

%% PLOT binary heatmap of nonzero stationary states by NODES

% ARGUMENTS pof function:
% fcn_plot_statsol_bin_hmap(stat_sol,prob_thresh,term_verts_inds_cell,nodes,sel_nodes,plot_param_settings,tight_subplot_flag,ranking_flag)
% stat_sol: vector of stationary solutions
% probability threshold for states to show (if left empty, all states shown)
prob_thresh=0.01;  % []; % 0.05;
% term_verts_cell: which subgraph to plot if there are disconnected ~
% nodes: name of nodes
% nodes to show. if none selected, then all nodes shown
sel_nodes=[]; % setdiff(2:numel(nodes)-1,[find(strcmp(nodes,{'Rb_b2'})) find(strcmp(nodes,{'p27_b2'}))]);
%
% plot_param_settings
% num_size_plot: font size of 0/1s on the heatmap
% hor_gap: horizontal gap between terminal SCCs, bottom_marg: bottom margin, left_marg: left margin
numsize_plot=22; fontsize=30; hor_gap=0.02; bottom_marg=0.31; left_marg=0.16; 
plot_param_settings=[numsize_plot fontsize hor_gap bottom_marg left_marg];
% want to use tight subplot? | order states by probability?
tight_subplot_flag='yes'; ranking_flag='yes';

% PLOT
figure('name','statsol_binary_heatmap')
statsol_binary_heatmap=fcn_plot_statsol_bin_hmap(stat_sol,prob_thresh,...
                        term_verts_cell, nodes,sel_nodes,plot_param_settings,tight_subplot_flag,ranking_flag);
% SAVE
% if <overwrite_flag> non-empty then existing file with same name is overwritten. 
overwrite_flag='yes'; 
% for <resolution> you can enter dpi value manually
% magnification=0.8; resolution_dpi=strcat('-r',num2str(magnification*get(0, 'ScreenPixelsPerInch')));
resolution_dpi='-r350'; fcn_save_fig('binary_heatmap_states',plot_save_folder,fig_file_type{3},overwrite_flag,resolution_dpi);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% STGs on subplots, with given parameter highlighted on each

% transitions in non-empty subgraph
popul_subgraphs=cellfun(@(x) sum(ismember(find(x0>0), x)), cell_subgraphs)>0;
subgraph_states=cell2mat(cell_subgraphs(popul_subgraphs)');
% subgraph_states=cell2mat(cell_subgraphs(~cellfun(@(x) isempty(x),term_verts_cell))');
% identify all transition rates that have corresponding transitions
par_inds_table = unique(stg_table(ismember(stg_table(:,1), subgraph_states) | ismember(stg_table(:,2), subgraph_states),3:4),'rows');
stg_table_subgraph = stg_table(ismember(stg_table(:,1), subgraph_states) | ismember(stg_table(:,2), subgraph_states),:);

for k=1:size(par_inds_table,1)
    param_freq(k)=sum(stg_table_subgraph(:,3)==par_inds_table(k,1) & stg_table_subgraph(:,4)==par_inds_table(k,2));
end
% top n most frequent transitions
[~,b]=maxk(param_freq,2);

% ARGUMENTS of function
% A_sparse,subgraph_index,term_verts_cell,cell_subgraphs,stat_sol,nodes: same as previous function
% selected_pars: parameters to highlight, either 'all', or numeric array [1 2 3]
% which node this parameter corresponds to?
sel_params=unique(par_inds_table(:,1)); 
% top n most frequent: sel_params=sel_params(b);
% all transition rates: unique(par_inds_table(:,1))'
% highest frequency: par_inds_table(b,1)'; % 
% up or down rate?
sel_params_up_down=arrayfun(@(x) par_inds_table(par_inds_table(:,1)==x,2)', sel_params,'un',0);

% stg_table: state transition table (generated by stg_table=fcn_build_stg_table(truth_table_filename,nodes))
plot_pars=[20 0.1 7 3 9]; % plot_pars=[fontsize,linewidth_val, arrowsize, default_markersize, highlight_markersize]
% parameters for highlighted transitions: color and width of highlighted edges
highlight_settings={'yellow',3}; 
% if using tight_subplots toolbox:
tight_subpl_flag='yes'; tight_subplot_pars=[0.06 0.02; 0.05 0.05; 0.05 0.05]; 

subgraph_index=find(cellfun(@(x) sum(ismember(find(x0>0), x)), cell_subgraphs)>0);

figure('name','STG select params');
plot_STG_sel_param(A_sparse,subgraph_index,term_verts_cell,cell_subgraphs,stat_sol,nodes,sel_params,sel_params_up_down,stg_table,...
    plot_pars,highlight_settings,tight_subpl_flag,tight_subplot_pars)

% SAVE
resolution_dpi='-r350'; fcn_save_fig('STG_subgraph4_highlighted_params_all',plot_save_folder,fig_file_type{1},overwrite_flag,resolution_dpi);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Parameter sensitivity analysis: one-dimensional parameter scans

% select the nodes whose parameters we want to scan in:
% all nodes that have actual transitions
% par_inds_table=unique(stg_table(:,[3 4]),'rows');
popul_subgraphs=cellfun(@(x) sum(ismember(find(x0>0), x)), cell_subgraphs)>0;
subgraph_states=cell2mat(cell_subgraphs(popul_subgraphs)');
% subgraph_states=cell2mat(cell_subgraphs(~cellfun(@(x) isempty(x),term_verts_cell))');
par_inds_table=unique(stg_table(ismember(stg_table(:,1), subgraph_states) | ismember(stg_table(:,2), subgraph_states),3:4),'rows');

% most common transitions
for k=1:size(par_inds_table,1)
    param_freq(k) = sum(stg_table(:,3)==par_inds_table(k,1) & stg_table(:,4)==par_inds_table(k,2));
end
% top n most frequent transitions
[~,top_freq_trans_rates]=maxk(param_freq,6);

% all rates that have corresponding transitions
scan_params=unique(par_inds_table(:,1))';
% all nodes that have transitions: unique(par_inds_table(:,1))'; 
% most frequent: par_inds_table(top_freq_trans_rates,1)';
% Zanudo: find(ismember(nodes,{'AKT','SGK1','TSC','FOXO3','BIM','BAD','mTORC1','PI3K','PRAS40'})); 
% all nodes that have transitions, except phenotypes: setdiff(unique(par_inds_table(:,1))',find(ismember(nodes,{'Apoptosis','Proliferation'})))
% selected nodes: find(ismember(nodes,{'AKT','SGK1','TSC','FOXO3','BIM','BAD','mTORC1'})); % 
scan_params_up_down=arrayfun(@(x) par_inds_table(par_inds_table(:,1)==x,2)', scan_params,'un',0); 
% how many transition rates we'll scan? sum(cellfun(@(x) numel(x),scan_params_up_down))
%
% num2cell(repelem([1 2],numel(scan_params),1),2)'; % both up and down rates
% num2cell(ones(1,numel(scan_params))); % only up 
% num2cell(2*ones(1,numel(scan_params))); % only down
% {[1 2], 1, [1 2]}; % manually selected

% top frequency trans rates:
% scan_params=par_inds_table(top_freq_trans_rates,1)'; 
% scan_params_up_down=arrayfun(@(x) par_inds_table( top_freq_trans_rates(par_inds_table(top_freq_trans_rates,1)==x),2)', scan_params,'un',0); 

% min and max of range of values; resolution of the scan; linear or logarithmic sampling
parscan_min_max = [1e-2 1e2]; n_steps=10; sampling_types={'log','linear'}; 

% FUNCTION for generating matrix of ordered values for the parameters to scan in
% [scan_par_table,scan_par_inds,~]= fcn_get_trans_rates_tbl_inds(scan_params,scan_params_up_down,transition_rates_table);
parscan_matrix=fcn_onedim_parscan_generate_matrix(scan_params,scan_params_up_down,nodes,sampling_types{1},parscan_min_max,n_steps);
% to set values in one/more columns manually (order is {3,1}->{3,2}->{4,1}->{4,2} etc)

% FUNCTION for 1-DIMENSIONAL PARSCAN
% ARGUMENTS
% stg table: generated above, state transition graph
% transition_rates_table:default values of transition rates
% initial conditions: x0
% stg_table: generated by stg_table=fcn_build_stg_table(truth_table_filename,nodes);
% [~,scan_par_inds,~]=fcn_get_trans_rates_tbl_inds(scan_params,scan_params_up_down,nodes);
tic;
[stationary_state_vals_onedimscan,stationary_node_vals_onedimscan,stationary_state_inds_scan]=...
    fcn_onedim_parscan_calc(stg_table,transition_rates_table,x0,nodes,parscan_matrix,scan_params,scan_params_up_down);
toc;

%% PLOT RESULTS of 1-by-1 parameter scan on heatmap/lineplot BY PARAMETERS

%%% FIRST PLOT TYPE: each subplot is a lineplot of node or state values as a function of a parameter's value, with a defined minimal variation

nonzero_states_inds=find(stat_sol>0);
% plot parameters
% [0.06 0.03],[0.03 0.03],[0.02 0.01]
height_width_gap=[0.08 0.03]; bott_top_marg =[0.05 0.05]; left_right_marg=[0.04 0.01];
% [fontsize_axes,fontsize_title,legend_fontsize,linewidth,params_tight_subplots(leave empty if not installed),model_name]
plot_param_settings={20,20,20,4,{height_width_gap bott_top_marg left_right_marg},model_name}; % plot_param_settings={12,14,[],model_name}; 
state_or_node_flags={'nodes','states'}; 
% cutoff for minimal variation to show a variable
diff_cutoff=0.1;
figure('name','onedim parscan by param')
[fig_filename,output_cell]=fcn_onedim_parscan_plot_by_params(state_or_node_flags{1},...
                                      stationary_node_vals_onedimscan,stationary_state_vals_onedimscan,...
                                      nonzero_states_inds,parscan_matrix,nodes,...
                                      scan_params,scan_params_up_down,... % selected parameters
                                      diff_cutoff,... % minimal variation for variable to be shown on plot
                                      plot_param_settings);
% SAVE figure
% resolution_dpi='-r350'; fcn_save_fig(strcat(fig_filename,'_r350'),plot_save_folder,fig_file_type{2},'overwrite',resolution_dpi);

%% PLOT RESULTS of 1-by-1 parameter scan on heatmap/lineplot BY VARIABLES

%%% SECOND PLOT TYPE: show stationary value/response coefficient of 1 variable (state or node) on 1 subplot, as a fcn of all relevant parameters
nonzero_states_inds=find(stat_sol>0);
sensit_cutoff=0.1; % minimal value for response coefficient (local sensitivity) or for the variation of node/state values
% nonzero states of the model
% nonzero_states=unique(cell2mat(stationary_state_inds_scan(:)'))';
% select parameters of plot
height_width_gap=[0.03 0.01]; bott_top_marg=[0.13 0]; left_right_marg=[0.04 0.04];
% [fontsize_axes,fontsize_title,params_tight_subplots(leave empty if not installed),model_name]
plot_param_settings={30,30,{height_width_gap bott_top_marg left_right_marg},model_name,'colorbar'}; 
% plot_param_settings={12,14,[],model_name}; 
% select type of plot
plot_types={{'lineplot','heatmap'} {'nodes','states'} {'values','sensitivity'}};
% if want to loop through all plot types: all_opts_perm=[[1 1 1]; unique([perms([1 1 2]); perms([2 2 1])],'rows'); [2 2 2]];
plot_type_options=[1 2 1];
figure('name',strjoin(arrayfun(@(x) plot_types{x}{plot_type_options(x)}, 1:numel(plot_type_options), 'un',0),'_'));
[resp_coeff,scan_params_sensit,scan_params_up_down_sensit,fig_filename]=fcn_onedim_parscan_plot_parsensit(plot_types,plot_type_options,...
                                                   stationary_node_vals_onedimscan,stationary_state_vals_onedimscan,...
                                                   nonzero_states_inds,parscan_matrix,nodes,...
                                                   scan_params,scan_params_up_down,... % transition rates to scan in
                                                   sensit_cutoff,plot_param_settings);

% <resp_coeffs> dimensions: (parameters, values,nodes), so eg. resp_coeffs(:,:,7)=resp. coeff values across param ranges for the 7th node

% SAVE figure
% resolution_dpi='-r350'; fcn_save_fig(strcat(fig_filename,'_cutoff',strrep(num2str(sensit_cutoff),'.','p')),plot_save_folder,fig_file_type{3},'overwrite',resolution_dpi);

%% multidimensional parameter scan: LATIN HYPERCUBE SAMPLING (random multidimensional sampling within given parameter ranges)

% WHICH PARAMETERS to scan simultaneously?
% all relevant parameters: 
% scan_params=unique(stg_table(:,3))'; scan_params_up_down=num2cell(repelem([1 2],numel(scan_params),1),2)'; 
% scan_params_up_down_sensit=arrayfun(@(x) scan_params_up_down_sensit{x}', 1:numel(scan_params_up_down_sensit),'un',0)

% PERFORM Latin Hypercube sampling (LHS) SAMPLING
sampling_types={'lognorm','linear','logunif'};
sampling_type=sampling_types{3};
% <lhs_scan_dim>: number of param sets
lhs_scan_dim=50;
% par_min_mean: minimum or in case of lognormal the mean of distribution. Can be a scalar or a vector, 
% if we want different values for different parameters
% max_stdev: maximum or in case of lognormal the mean of distribution. 
% Can be a scalar or a vector, if we want different values for different parameters
%
% for 'lognorm' and 'logunif' provide the LOG10 value of desired mean/min and stdev/max!!, ie. -2 means a mean of 0.01
par_min_mean=-2; % repmat(1.5,1,numel(cell2mat(scan_params_up_down(:)))); par_min_mean(4)=3; 
max_stdev=2; % repmat(0.5,1,numel(cell2mat(scan_params_up_down(:))));
[all_par_vals_lhs,stat_sol_nodes_lhs_parscan,stat_sol_states_lhs_parscan,stat_sol_states_lhs_parscan_cell]=... % outputs
    fcn_multidim_parscan_latinhypcube(par_min_mean,max_stdev,sampling_type,lhs_scan_dim, ...
                                            scan_params_sensit,scan_params_up_down_sensit, ... % transition rates
                                            transition_rates_table,stg_table,x0,nodes);

% OUTPUTS:
% <all_par_vals_lhs>: table of parameter sets
% <stat_sol_nodes_lhs_parscan>: stationary values of nodes
% <stat_sol_states_lhs_parscan>: stationary values of states
% <stat_sol_states_lhs_parscan_cell>: indices of non-zero states, if they are always the same the variable is a 1-row array, 
% otherwise a cell with the non-empty elements stat_sol_states_lhs_parscan_cell{k,1} contain values of nonzero
% states, elements stat_sol_states_lhs_parscan_cell{k,2} their indices 

%% SCATTERPLOTS of STATE or NODE values as a function of the selected parameters, with the trendline shown (average value per parameter bin)

% sel_nodes=[6 10 11 12 13 14 15];
for var_ind=scan_params_sensit
% find(strcmp(nodes,'CHEK1')); % which STATE or NODE to plot
% <all_par_vals_lhs>: parameter sets
% [number_bins_for_mean,trendline_width,axes_fontsize,index nonzero states]
param_settings = [50 4 24 size(stat_sol_states_lhs_parscan_cell)]; 
% STATES or NODES? <scan_values>: values to be plotted
scan_values=stat_sol_nodes_lhs_parscan; % stat_sol_states_lhs_parscan
% PLOT
sampling_type=sampling_types{3}; % sampling_types={'lognorm','linear','logunif'};
file_name_prefix=strcat('LHS_parscan_trend_',nodes{var_ind}); 
figure('name',nodes{var_ind})
fcn_multidim_parscan_scatterplot(var_ind,all_par_vals_lhs,scan_values,...
        scan_params_sensit,scan_params_up_down_sensit,nodes,sampling_type,param_settings)
end

resolution_dpi='-r350'; fcn_save_fig(file_name_prefix,plot_save_folder,fig_file_type{3},'overwrite',resolution_dpi);

%% calculating & plotting (heatmap) correlations between variables 

% ARGUMENTS
% stat_sol_nodes_lhs_parscan: values from parameter sampling
% nodes: name of ALL nodes
% sel_nodes: name of selected nodes (pls provide in ascending order)
% fontsize: ~ for labels and titles (displaying correlation)
% HEATMAPS of correlations between selected variables
sel_nodes=[]; % 3:numel(nodes); % scan_params_sensit
plot_settings=[10 8 8]; % [fontsize on plot, fontsize on axes/labels]
plot_type_flag={'var_var','heatmap'}; % this is plotting the heatmap of correlations between variables
figure('name',strjoin(plot_type_flag))
[varvar_corr_matr,p_matrix_vars]=fcn_multidim_parscan_parvarcorrs(plot_type_flag,all_par_vals_lhs,stat_sol_nodes_lhs_parscan,...
                                            nodes,sel_nodes,[],[],[],plot_settings);

% resolution_dpi='-r350'; fcn_save_fig(strcat(strjoin(plot_type_flag,'_'),'_corrs'),plot_save_folder,fig_file_type{1},'overwrite',resolution_dpi);
                   
%% scatterplots of selected variables [i,j]: var_i VS var_j

sel_nodes=scan_params_sensit; plot_settings=[10 12 12]; % [fontsize_axes, fontsize_titles]
plot_type_flag={'var_var','scatter'}; % this is plotting the scatterplots of variables with correlation values
figure('name',strjoin(plot_type_flag))
fcn_multidim_parscan_parvarcorrs(plot_type_flag,all_par_vals_lhs,stat_sol_nodes_lhs_parscan,...
                                    nodes,sel_nodes,[],[],[],plot_settings);

resolution_dpi='-r350'; 
fcn_save_fig(strcat(strjoin(plot_type_flag,'_'),'_scatterplot'),plot_save_folder,fig_file_type{3},'overwrite',resolution_dpi);

%% linear or lin-log regression of VARIABLES as fcn of PARAMETERS: VARIABLE=f(PARAMETER), the function plots R squared

plot_type_flag={'par_var','heatmap','r_sq'}; % {'par_var','heatmap'/'lineplot','r_sq'/'slope'}
sel_nodes=numel(nodes); % 3:numel(nodes); % scan_params_sensit;
% plot_settings=[fontsize,maximum value for heatmap colors], if plot_settings(2)=NaN, then max color automatically selected
plot_settings=[20 20 0.5]; 
% if regression type is 'linlog', then the fit is y = a + b*log10(x)
regr_types={'log','linear'}; % linlog recommended if parameter values log-uniformly distributed in sampling
figure('name',strjoin(plot_type_flag))
[r_squared,slope_intercept]=fcn_multidim_parscan_parvarcorrs(plot_type_flag,all_par_vals_lhs,stat_sol_nodes_lhs_parscan,...
                                 nodes,sel_nodes,... % which nodes
                                 scan_params_sensit,scan_params_up_down_sensit, ... % parameters (CAREFUL that they are the same as in LHS!)
                                 regr_types{1},plot_settings);

resolution_dpi='-r350'; fcn_save_fig(strjoin(plot_type_flag,'_'),plot_save_folder,fig_file_type{3},'overwrite',resolution_dpi)

%% Quantify importance of parameters from LHS by a regression tree

% predictor importance values: look into arguments of <fitrtree> in MATLAB documentation to modulate regression tree
% for STATES or NODES?
scan_values=stat_sol_nodes_lhs_parscan; % stat_sol_states_lhs_parscan
sel_nodes=scan_params_sensit; % STATES or NODES to be analyzed
% names of selected transition rates and their predictor importance values
% [~,~,predictor_names] = fcn_get_trans_rates_tbl_inds(scan_params,scan_params_up_down,nodes);
% predictorImportance_vals=cell2mat(arrayfun(@(x) predictorImportance(...
%     fitrtree(all_par_vals_lhs,scan_values(:,x),'PredictorNames',predictor_names)), sel_nodes,'un',0)');

% CALCULATE and PLOT predictor importance
plot_type_flags={'line','bar'};
figure('name','regression_tree_pred_import')
% [predictor_names,predictorImportance_vals]
[~,predictorImportance_vals] = fcn_multidim_parscan_predictorimport(scan_params_sensit,scan_params_up_down_sensit,...
                                                all_par_vals_lhs,scan_values,nodes,sel_nodes,plot_type_flags{2});
                                            
% resolution_dpi='-r350'; fcn_save_fig('regression_tree_pred_import',plot_save_folder,fig_file_type{3},'overwrite',resolution_dpi)
   
%% Sobol total sensitivity metric                                            

% On Sobol total sensitivity index see: https://en.wikipedia.org/wiki/Variance-based_sensitivity_analysis
% This metric indicates how much of the total variance in a variable is due to variation in a given parameter
% We calculate here the usual numerical approximation of analytical equivalent from Monte Carlo sampling.
% From the LHS sampling above we take the matrices of parameter sets and variable values:
% [parameter sets, variable values]: [all_par_vals_lhs,stat_sol_nodes_lhs_parscan]

% select only parameters with an R-squared over some value
% careful to use same parameters as for R^2 calculations!!!
[par_ind_table,sequential_indices_lhs,~] = fcn_get_trans_rates_tbl_inds(scan_params_sensit,scan_params_up_down_sensit,nodes);
% R^2>=0.1
r_sq_thresh=0.05;
par_ind_table_filtered=par_ind_table(sum(r_squared>r_sq_thresh)>0,:);
scan_params_filtered=unique(par_ind_table_filtered(:,1))'; % scan_params_unique=unique(par_ind_table(sum(r_squared>0.2)>0,1))'; 
scan_params_up_down_filtered=arrayfun(@(x) par_ind_table_filtered(par_ind_table_filtered(:,1)==x,2)', scan_params_filtered,'un',0);

% Sobol total sensitivity: calculated for one variable at a time
% selected nodes to display
sel_nodes=[]; % 3:numel(nodes); % scan_params_sensit;
% setdiff(1:numel(nodes),find(sum(cell2mat(arrayfun(@(x) strcmp(nodes,x), {'cc','KRAS','CDC25B'},'un',0)'))));
% sel_vars=[]; % if left empty, all nodes/states are analyzed
sample_size=[]; % if left empty, the sample size is half of the original param scan <all_par_vals_lhs>
% PLOT SETTINGS: [fontsize_plot,fontsize_axes,fontsize_title, min_color(optional), max_color(opt), progress_calcul_every_x_% (opt)];
plot_settings=[20 30 30 NaN NaN 10];
var_types={'node','state'}; % analysis for states or nodes
% to calculate Sobol total sensitivity we need <sample_size*numel(scan_params_up_down)> evaluations of the model
figure('name','sobol sensitivity index')
sobol_sensit_index=fcn_multidim_parscan_sobol_sensit_index([],var_types{1},...
                      all_par_vals_lhs,stat_sol_nodes_lhs_parscan,stat_sol_states_lhs_parscan,...
                      sample_size,... % # of calculations per parameter
                      sequential_indices_lhs,... % this is indices of transition rates in the original LHS
                      scan_params_filtered,scan_params_up_down_filtered,...% scan_params_sensit,scan_params_up_down_sensit
                      stg_table,transition_rates_table,x0,nodes,sel_nodes,plot_settings);

% if already calculated <sobol_sensit_index> and only want to plot results, provide <sobol_sensit_index> as FIRST argument 
fcn_multidim_parscan_sobol_sensit_index(sobol_sensit_index,var_types{1},[],[],[],[],...
                       sequential_indices_lhs,scan_params_filtered,scan_params_up_down_filtered,[],[],[],nodes,sel_nodes,plot_settings);

% SAVE
% magnification=0.8; resolution_dpi=strcat('-r',num2str(magnification*get(0,'ScreenPixelsPerInch')));
resolution_dpi='-r350'; fcn_save_fig('sobol_sensitivity_index',plot_save_folder,fig_file_type{3},'overwrite',resolution_dpi)

%% PARAMETER FITTING

% define parameters to vary (predictor_names)
% sensitive parameters identified by 1-dimensional param scan: scan_params_sensit,scan_params_up_down_sensit
[~,~,predictor_names]=fcn_get_trans_rates_tbl_inds(scan_params_sensit,scan_params_up_down_sensit,nodes); 
% define data vector (generate some data OR load from elsewhere)
sel_param_vals=lognrnd(1,1,1,numel(predictor_names)); % abs(normrnd(1,0.5,1,numel(predictor_names)));
transition_rates_table=fcn_trans_rates_table(nodes,'uniform',[],[],predictor_names,sel_param_vals);
y_data=fcn_calc_init_stat_nodevals(x0,split_calc_inverse(fcn_build_trans_matr(stg_table,transition_rates_table,''),stg_sorting_cell,...
                                   transition_rates_table,x0),'x0');

% create functions that calculate sum of squared deviations & values of
% variables (composed of different fcns) - RERUN THIS if you want to fit to new data or new non-fitted transition rates!!
[fcn_statsol_sum_sq_dev,fcn_statsol_values]=fcn_simul_anneal(y_data,x0,stg_table,stg_sorting_cell,nodes,predictor_names);

% FITTING by simulated annealing (look at arguments in anneal/anneal.m)
% initial guess for parameters
init_vals=lognrnd(0,2,size(predictor_names)); init_error=fcn_statsol_sum_sq_dev(init_vals);
% initial value of model nodes
y_init=fcn_calc_init_stat_nodevals(x0,...
    split_calc_inverse(fcn_build_trans_matr(stg_table,fcn_trans_rates_table(nodes,'uniform',[],[],predictor_names,init_vals),''),stg_sorting_cell,...
    transition_rates_table,x0),'');

% simulated annealing with existing algorithm anneal/anneal.m, with modifications in script: 
% 1) defined counter=0 before while loop
% 2) inserted <T_loss(counter,:)=[T oldenergy];> at line 175, defined <T_loss> as 3rd output
% arguments for algorithm need to be provided as a structure: 
% struct('Verbosity',2, 'StopVal', 0.01, 'StopTemp',1e-8) % stopping error value, stopping temperature parameter
fitting_arguments=struct('Verbosity',2, 'StopVal', 0.002);
tic; [optim_par_vals,best_error,T_loss]=anneal(fcn_statsol_sum_sq_dev,init_vals,fitting_arguments); toc 
% 20-40 mins for 15var KRAS model

% output with fitted parameters
y_optim_param=fcn_calc_init_stat_nodevals(x0,...
    split_calc_inverse(fcn_build_trans_matr(stg_table,fcn_trans_rates_table(nodes,'uniform',[],[],predictor_names,optim_par_vals),''),stg_sorting_cell,...
    transition_rates_table,x0),'');
% plot data, initial values, optimized values
data_init_optim=[y_data; y_init; y_optim_param]; min_val=min(min(data_init_optim(:,3:end))); max_val=max(max(data_init_optim(:,3:end)));

figure('name','param fitting'); sel_nodes=3:numel(nodes);
% PLOT fitting process
thres_ind=size(T_loss,1); % thres_ind=find(T_loss(:,2)<1e-2,1); 
vars_show=2; % 1=temperature, 2=error
% figure('name','simul anneal')
fig_subpl1=subplot(1,2,1);
plot(1:thres_ind, T_loss(1:thres_ind,vars_show),'LineWidth',4); xlim([0 thres_ind]); if init_error/best_error>30; set(gca,'yscale','log'); end
label_ticks_fontsize=24; label_fontsize=30; set(gca,'FontSize',label_ticks_fontsize);
legend_strs={'temperature', 'sum of squared error'}; ylabel(legend_strs{vars_show},'FontSize',label_fontsize); grid on; 
xlabel('number of iterations','FontSize',label_fontsize); % title('Parameter fitting by simulated annealing','FontSize',22)
% 
fig_subpl2=subplot(1,2,2); 
barplot_gca=barh(data_init_optim(:,sel_nodes)'); set(fig_subpl2,'ytick',1:numel(sel_nodes)); 
legend({'data','initial guess','optimized'},'FontSize',label_fontsize)
set(gca,'FontSize',label_ticks_fontsize); 
xticklabels=get(gca,'xtick'); set(fig_subpl2,'xticklabel',xticklabels,'FontSize',label_ticks_fontsize);
xlabel('stationary probabilities','FontSize',label_fontsize); 
set(fig_subpl2,'yticklabel','');  % strrep(nodes(sel_nodes),'_','\_'),'FontSize',label_fontsize

resolution_dpi='-r350';
fcn_save_fig(strcat('simulated_annealing_',num2str(numel(predictor_names)),'fittingpars'),plot_save_folder,fig_file_type{3},'overwrite',resolution_dpi)

% mean absolute error: mean(abs(y_data - fcn_statsol_values(optim_par_vals)))
% distance of fitted params from true values: abs(optim_par_vals - sel_param_vals)./sel_param_vals