% This file contains the commands to run the functions calculating the
% stationary solution of stoch logical models, plot results and perform parametric analysis

% go to the folder of the file
editor_service = com.mathworks.mlservices.MLEditorServices; editor_app = editor_service.getEditorApplication;
active_editor = editor_app.getActiveEditor; storage_location = active_editor.getStorageLocation;
file = char(storage_location.getFile); path_to_toolbox = fileparts(file); cd(path_to_toolbox);

%% ADD FUNCTIONS to PATH

addpath('functions/') 
% using these publicly available MATLAB toolboxes:
% needed for PLOTS:
addpath('heatmaps') % https://mathworks.com/matlabcentral/fileexchange/24253-customizable-heat-maps
addpath('redblue'); % https://mathworks.com/matlabcentral/fileexchange/25536-red-blue-colormap (for redblue colormaps)

% optional for generating and saving plots
% for subplots with smaller gaps:
addpath('tight_subplot') % https://mathworks.com/matlabcentral/fileexchange/27991-tight_subplot-nh-nw-gap-marg_h-marg_w

% export figures as EPS or PDF as they appear on plots:
% Optional. https://mathworks.com/matlabcentral/fileexchange/23629-export_fig
export_fig_name=dir('altman*'); export_fig_name=export_fig_name.name; addpath(genpath(export_fig_name)) 

% optional for PARAMETER FITTING by simulated annealing:
addpath('anneal') % https://mathworks.com/matlabcentral/fileexchange/10548-general-simulated-annealing-algorithm

% optional defaults settings: title font weight normal, docked figures
set(0,'DefaultAxesTitleFontWeight','normal'); set(0,'DefaultFigureWindowStyle','docked');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% model set-up

% models can be defined 
% A) by entering the list of nodes and their
% corresponding rules as a cell of strings, using MATLAB logical notation ('&', '|', '~', '(', ')'),
% for instance:
% nodes = {'cc','kras', 'dna_dam', 'chek1', 'mk2', 'atm_atr', 'hr','cdc25b', 'g2m_trans', 'cell_death'};
% 
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
model_name='mammalian_cc'; % kras15vars

% where to save figures
save_folder=strcat('doc/sample_plots/',model_name,'/');

% OR B) the model can be read in from an existing BOOLNET file
[nodes,rules]=fcn_bnet_readin('model_files/traynard2016_mammalian_cellcycle.bnet'); % krasmodel10vars.bnet

% once we have the list of nodes and their logical rules, we can check if
% all variables referred to by rules are found in the list of nodes:
fcn_nodes_rules_cmp(nodes,rules)

% if yes, we generate a function file, which will create the model
truth_table_filename='fcn_truthtable.m'; fcn_write_logicrules(nodes,rules,truth_table_filename)

% from the model we generate the STG table, that is independent of values of transition rates
tic; stg_table=fcn_build_stg_table(truth_table_filename,nodes); toc

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

% density of transition matrix A
% nnz(A_sparse)/numel(A_sparse)

%% calculate steady states with kernel/nullspace calculation functions

% defining an initial condition
n_nodes=numel(nodes); truth_table_inputs=rem(floor([0:((2^n_nodes)-1)].'*pow2(0:-1:-n_nodes+1)),2);
% define some nodes with a fixed value and a probability <dom_prob>: states
% that satisfy this condition will have a total initial probability of
% <dom_prob>, the other states 1-dom_prob

initial_fixed_nodes={'CycE','CycA','CycB','Cdh1','Rb_b1','Rb_b2','p27_b1','p27_b2'}; initial_fixed_nodes_vals=[0 0 0 1 1 1 1 1];
% KRAS model, WT: {'cc','KRAS'}, [1 0]. Mutant: {'cc','KRAS'}, [1 0]
% for cell cycle model initial state: CycE=0 & CycA=0 & CycB=0 & Cdh1=1 & % Rb=1 & p27=1, meaning 
% {CycE CycA CycB Cdh1 Rb_b1 Rb_b2 p27_b1 p27_b2} [0 0 0 1 1 1 1 1]
% what is the probability of this state, (eg. dom_prob=0.8, ie. 80% probability)
dom_prob=1;
% if <random> the probability is randomly distributed among states, if <uniform> uniformly
distrib_types={'random','uniform'}; 
% if plot_flag non-empty, we get a bar plot of initial values
plot_flag='';
% function assigns a probability of <dom_prob> to the states with the fixed nodes having the defined values
x0=fcn_define_initial_states(initial_fixed_nodes,initial_fixed_nodes_vals,dom_prob,nodes,distrib_types{1},plot_flag);

% completely random initial condition: 
% x0=zeros(2^n_nodes,1); x0=rand(1,size(truth_table_inputs,1))'; x0=x0/sum(x0);
% completely uniform initial condition
% x0=ones(1,2^n_nodes)/2^n_nodes;

%% CALCULATE STATIONARY STATE
% ARGUMENTS:
% transition matrix: A
% table of transition rates: transition_rates_table
% initial conditions: x0
tic; [stat_sol,term_verts_cell,cell_subgraphs]=split_calc_inverse(A_sparse,transition_rates_table,x0); toc
% OUTPUTS
% stat_sol: stationary solution for all the states
% term_verts_cell: index of nonzero states. If the STG is disconnected the nonzero states corresp to these disconn subgraphs are in separate cells
% cell_subgraphs: indices of states belonging to disconnected subgraphs (if any)
% 
% probabilities by subgraph:
% arrayfun(@(x) sum(stat_sol(cell2mat(term_verts_cell{x}))), 1:numel(term_verts_cell))

% nonzero states can be quickly queried by:
stat_sol(stat_sol>0)' % probability values of nonzero states
truth_table_inputs(stat_sol>0,:) % logical states that are nonzero

% sum the probabilities of nonzero states by nodes, both for the initial condition and the stationary solution
% ARGUMENTS
% initial conditions: x0
% % stat_sol: stationary solution for all the states
% nodes: list of nodes
[stationary_node_vals,init_node_vals]=fcn_calc_init_stat_nodevals(x0,stat_sol);

% can check results by doing direct matrix exponentiation, for systems up to ~10 nodes, but note this will not be stationary if there are cycles.
% x_sol=((x0')*A_sparse^1e5)';
% stationary_node_vals=x_sol'*truth_table_inputs;
% 
% checked with MaBoSS simuls, results are identical (up to 1% dev.) as they have to be
% comparing with simulation of mammalian cell cycle model with 12
% nodes: look in folder <sample_plots/mammalian_cc/maboss_files_plots>

%% PLOTTING RESULTS

% PLOT A/K and stat solutions
% ARGUMENTS
% matrix_input: [], K_sparse or A_sparse (kinetic or transition matrix)
% min_max_col: minimum and max color for heatmap
% fontsize: [font size on the heatmap, title font size for stationary solutions]
% barwidth_states_val: width of the bars for bar plot of stationary solutions of states
% sel_nodes: nodes to show. If left empty, all nodes are shown
% nonzero_flag: minimal value for probability to display - if this is non-empty, only plot nonzero states, useful for visibility if there are many states
sel_nodes=[]; min_max_col=[0 1]; barwidth_states_val=0.8;fontsize=[10 20]; % fontsize_hm,fontsize_stat_sol
plot_settings = [fontsize barwidth_states_val min_max_col]; prob_thresh=0.01;
% WARNING!!! if more than 12 nodes, generating the figure for A/K can be time-consuming
matrix_input=A_sparse;
figure('name','A_K_stat_sol')
fcn_plot_A_K_stat_sol(matrix_input, nodes, sel_nodes, stat_sol, x0, plot_settings,prob_thresh)

% SAVE
% enter any string for the last argument to overwrite existing plot!!
if exist(save_folder,'dir')==0; mkdir(save_folder); end
fig_file_type={'.png','.eps','.pdf','.jpg','.tif'}; if ~isempty(matrix_input); matrix_input_str='_with_matrix'; else; matrix_input_str=''; end
% if <overwrite_flag> non-empty then existing file with same name is overwritten. 
overwrite_flag='yes'; 
% for <resolution> you can enter dpi value manually, if left empty then it is manually set
magnification=0.8; resolution_dpi=strcat('-r',num2str(magnification*get(0, 'ScreenPixelsPerInch')));
fcn_save_fig(strcat('single_solution_states_nodes_stat_sol',matrix_input_str),save_folder,fig_file_type{2},overwrite_flag,resolution_dpi);

%% PLOT stationary solutions (without A/K matrix)

% nonzero_flag: if non-empty, only the nonzero states with a probability above this value are shown
sel_nodes=[]; prob_thresh=0.01; barwidth_states_val=0.8; % for 3 nonzero states ~0.8 is a good value
fontsize=[9 20]; % [fontsize_y_axis_states,fontsize_x_axes_and_titles]
plot_settings=[fontsize barwidth_states_val]; matrix_input=[];
figure('name','stat_sol')
fcn_plot_A_K_stat_sol(matrix_input, nodes, sel_nodes, stat_sol, x0, plot_settings,prob_thresh)

% SAVE
% enter any string for the last argument to overwrite existing plot!!
% if <overwrite_flag> non-empty then existing file with same name is overwritten. 
overwrite_flag='yes'; 
% for <resolution> you can enter dpi value manually, if left empty then it is manually set
magnification=0.8; resolution_dpi=strcat('-r',num2str(magnification*get(0, 'ScreenPixelsPerInch')));
fcn_save_fig(strcat('single_solution_states_nodes_stat_sol'),save_folder,fig_file_type{2},overwrite_flag,resolution_dpi);

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
numsize_plot=12; fontsize=12; hor_gap=0.01; bottom_marg=0.15; left_marg=0.04; 
plot_param_settings=[numsize_plot fontsize hor_gap bottom_marg left_marg];
% want to use tight subplot? | order states by probability?
tight_subplot_flag='yes'; ranking_flag='yes';
% PLOT
figure('name','statsol_binary_heatmap')
% inputting terminal vertices if there are multiple subgraphs and in some
% of them there are fixed points, in others a cyclic attractor: 
% for mammalian cell cycle model: [term_verts_cell{1} term_verts_cell{2}]
statsol_binary_heatmap=fcn_plot_statsol_bin_hmap(stat_sol,prob_thresh,...
                            [term_verts_cell{1} term_verts_cell{2}],... % if providing a single cell: term_verts_cell{nonempty_subgraph}
                            nodes,sel_nodes,plot_param_settings,tight_subplot_flag,ranking_flag);
% SAVE
% if <overwrite_flag> non-empty then existing file with same name is overwritten. 
overwrite_flag='yes'; 
% for <resolution> you can enter dpi value manually, if left empty then it is manually set
magnification=0.8; resolution_dpi=strcat('-r',num2str(magnification*get(0, 'ScreenPixelsPerInch')));
fcn_save_fig('binary_heatmap_states',save_folder,fig_file_type{2},overwrite_flag,resolution_dpi);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plot of STG (state transition graph) with source and sink (terminal) vertices highlighted

% NOTE: for models larger than ~10-12 nodes generating these plots can be very time consuming!

% calculate the sum of probabilities per disconnected subgraphs
probs_by_subgraph=arrayfun(@(x) sum(stat_sol(cell2mat(term_verts_cell{x}))), 1:numel(term_verts_cell)); % [~,subgraph_index]=max(probs_by_subgraph); 

% ARGUMENTS of function: 
% plot_STG(A_sparse,subgraph_index,term_verts_cell,cell_subgraphs,stat_sol,plot_settings,title_str,source_color)
%
% A_sparse: transition matrix for entire STG or selected subgraph
% subgraph_index: index of the subgraph to plot. If empty, entire STG is plotted
% term_verts_cell: cell of terminal vertices 
% cell_subgraphs: cell of vertices belonging to different subgraphs
% stat_sol: stationary values of probabilities of different states (vertices of STG)
% plot settings
% the entries are: [fontsize,linewidth_val,arrowsize,default_markersize for all vertices,marker size for terminal and sink vertices]
% size of source and sink states is proportional to their probability values
plot_settings=[20 1 7 3 9]; 
% color of vertices that have only outgoing edges
source_color='green';
figure('name','subgraph')
subgraph_index=1; title_str=strcat('subgraph #', num2str(subgraph_index));
subplot(1,2,1); stg_plot=plot_STG(A_sparse,subgraph_index,term_verts_cell,cell_subgraphs,stat_sol,plot_settings,title_str,source_color);
subgraph_index=2; title_str=strcat('subgraph #', num2str(subgraph_index)); 
subplot(1,2,2); plot_STG(A_sparse,2,term_verts_cell,cell_subgraphs,stat_sol,plot_settings,title_str,source_color)
% SAVE
magnification=0.8; resolution_dpi=strcat('-r',num2str(magnification*get(0, 'ScreenPixelsPerInch')));
fcn_save_fig('STG_subgraph',save_folder,fig_file_type{1},overwrite_flag,resolution_dpi);

%% STGs on subplots, with given parameter highlighted on each

% if STG table doesn't exist,  generate with <[stg_table,~,~]=fcn_build_stg_table(truth_table_filename,nodes);>

selected_pars=[1 3 4 5 6 11 9 8]; % parameters to highlight, either 'all', or numeric array [1 2 3]
plot_pars=[20 0.1 7 3 6]; % plot_pars=[fontsize,linewidth_val, arrowsize, default_markersize, highlight_markersize]
% parameters for highlighted transitions: color and width of corresp edges
highlight_settings={'yellow',3}; 
% if using tight_subplots toolbox:
tight_subpl_flag='yes'; tight_subplot_pars=[0.06 0.02; 0.05 0.05; 0.05 0.05]; 
% cropping plot (optional, for better visibility)
limits=[-4 5;-5 5]; 

figure('name','STG select params');
plot_STG_sel_param(A_sparse,counter,nodes,cell_subgraphs,selected_pars,stg_table,...
    plot_pars,highlight_settings,limits,tight_subpl_flag,tight_subplot_pars)
% show all params that have an effect
% plot_STG_sel_param(A_sparse,counter,nodes,cell_subgraphs,'all',stg_table,plot_pars,highlight_settings,'',tight_subpl_flag,tight_subplot_pars)

%% only one graph: in this case this is a subgraph of the entire STG, but can be the entire STG too

subgraph_index=2;
sample_size=1e3; sample_nodes=unique([cell2mat(term_verts_cell{subgraph_index})', datasample(cell_subgraphs{subgraph_index},sample_size)]); 
A_sub=A_sparse(sample_nodes,sample_nodes); 

% B=A_sparse(cell_subgraphs{counter},cell_subgraphs{counter});
% transition rates that occur in given subgraph
A_state_transitions_inds=fcn_stg_table_subgraph(stg_table,cell_subgraphs,subgraph_index);
par_inds_table=unique(A_state_transitions_inds(:,[3 4]),'row'); 
scan_params=unique(par_inds_table(:,1))'; scan_params_up_down=arrayfun(@(x) par_inds_table(par_inds_table(:,1)==x,2)', scan_params,'un',0); 
% sequential index of transition rates
[~,param_inds_sequential,~]=fcn_get_trans_rates_tbl_inds(scan_params,scan_params_up_down,nodes); 

% 2nd, 4th arguments left empty
figure('name','STG select params subgraph');
plot_STG_sel_param(A_sub,'',nodes,'',datasample(param_inds_sequential,5,'replace',false),...
    A_state_transitions_inds,default_settings,highlight_settings,[],'yes',tight_subplot_pars)

fcn_save_fig('STG_subgraph_highlighted_params',save_folder,fig_file_type{1},'');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
% stg_table: generated by stg_table=fcn_build_stg_table(truth_table_filename,nodes);
% [~,scan_par_inds,~]=fcn_get_trans_rates_tbl_inds(scan_params,scan_params_up_down,nodes);
[stationary_state_vals_onedimscan,stationary_node_vals_onedimscan,stationary_state_inds_scan]=...
    fcn_onedim_parscan_calc(stg_table,transition_rates_table,x0,nodes,parscan_matrix,scan_params,scan_params_up_down);

%% PLOT RESULTS of 1-by-1 parameter scan on heatmap/lineplot

%%% SECOND PLOT TYPE: show the stationary value or response coefficient of 1 variable or state on 1 subplot, as a fcn of all relevant parameters
nonzero_states_inds=find(stat_sol>0);
sensit_cutoff=0.1; % minimal value for response coefficient (local sensitivity)
% nonzero states of the model
% nonzero_states=unique(cell2mat(stationary_state_inds_scan(:)'))';
% select parameters of plot
height_width_gap=[0.03 0.03]; bott_top_marg =[0.13 0.1]; left_right_marg=[0.05 0.05];
% [fontsize_axes,fontsize_title,params_tight_subplots(leave empty if not installed),model_name]
plot_param_settings={12,14,{height_width_gap bott_top_marg left_right_marg},model_name}; % plot_param_settings={12,14,[],model_name}; 
% select type of plot
plot_types={{'lineplot','heatmap'} {'nodes','states'} {'values','sensitivity'}};
% if want to loop through all ploty types: all_opts_perm=[[1 1 1]; unique([perms([1 1 2]); perms([2 2 1])],'rows'); [2 2 2]];
plot_type_options=[2 2 2];
figure('name',strjoin(arrayfun(@(x) plot_types{x}{plot_type_options(x)}, 1:numel(plot_type_options), 'un',0),'_'));
[resp_coeff,scan_params_sensit,scan_params_up_down_sensit,fig_filename]=fcn_onedim_parscan_plot_parsensit(plot_types,plot_type_options,...
                                                           stationary_node_vals_onedimscan,stationary_state_vals_onedimscan,...
                                                           nonzero_states_inds,parscan_matrix,nodes,...
                                                           scan_params,scan_params_up_down,... % transition rates to scan in
                                                           sensit_cutoff,plot_param_settings);

% <resp_coeffs> dimensions: (parameters, values,nodes), so eg. resp_coeffs(:,:,7) are the resp. coeff values across the param ranges for the 7th node

% SAVE figure
fcn_save_fig(strcat(fig_filename),save_folder,fig_file_type{1},'overwr');

%% multidimensional parameter scan: LATIN HYPERCUBE SAMPLING (random multidimensional sampling within given parameter ranges)

% WHICH PARAMETERS to scan simultaneously?
% all relevant parameters: 
% scan_params=unique(stg_table(:,3))'; scan_params_up_down=num2cell(repelem([1 2],numel(scan_params),1),2)'; 
% scan_params_up_down_sensit=arrayfun(@(x) scan_params_up_down_sensit{x}', 1:numel(scan_params_up_down_sensit),'un',0)

% PERFORM Latin Hypercube sampling (LHS) SAMPLING
sampling_types={'lognorm','linear','logunif'};
sampling_type=sampling_types{3};
% <lhs_scan_dim>: number of param sets
lhs_scan_dim=1e3;
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

sel_nodes=[6 10 11 12 13 14 15];
for k=sel_nodes
var_ind=k; % find(strcmp(nodes,'CHEK1')); % which STATE or NODE to plot
% <all_par_vals_lhs>: parameter sets
param_settings = [50 4 16 stat_sol_states_lhs_parscan_cell]; % [number_bins_for_mean,trendline_width,axes_fontsize,index nonzero states]
% STATES or NODES? <scan_values>: values to be plotted
scan_values=stat_sol_nodes_lhs_parscan; % stat_sol_states_lhs_parscan
% PLOT
sampling_type=sampling_types{3}; % sampling_types={'lognorm','linear','logunif'};
file_name_prefix=strcat('LHS_parscan_trend_',nodes{var_ind}); figure('name',file_name_prefix)
fcn_multidim_parscan_scatterplot(var_ind,all_par_vals_lhs,scan_values,...
        scan_params_sensit,scan_params_up_down_sensit,...
        nodes,sampling_type,param_settings)
end

% full_filename_with_path=fcn_save_fig(file_name_prefix,save_folder,fig_file_type{1},'y');

%% calculating & plotting (heatmap) correlations between variables OR between params and variables by linear/logist regression

% ARGUMENTS
% stat_sol_nodes_lhs_parscan: values from parameter sampling
% nodes: name of ALL nodes
% sel_nodes: name of selected nodes (pls provide in ascending order)
% fontsize: ~ for labels and titles (displaying correlation)
% HEATMAPS of correlations between selected variables
sel_nodes=3:15; plot_settings=[15 16]; % [fontsize on plot, fontsize on axes/labels]
plot_type_flag={'var_var','heatmap'}; % this is plotting the heatmap of correlations between variables
figure('name',strjoin(plot_type_flag))
[varvar_corr_matr,p_matrix_vars]=fcn_multidim_parscan_parvarcorrs(plot_type_flag,all_par_vals_lhs,stat_sol_nodes_lhs_parscan,...
                                            nodes,sel_nodes,[],[],[],plot_settings);

full_filename_with_path=fcn_save_fig(strjoin(plot_type_flag,'_'),save_folder,fig_file_type{1});
                   
%% scatterplots of selected variables [i,j]: var_i VS var_j

sel_nodes=10:15; plot_settings=[10 12]; % [fontsize_axes, fontsize_titles]
plot_type_flag={'var_var','scatter'}; % this is plotting the scatterplots of variables with correlation values
figure('name',strjoin(plot_type_flag))
fcn_multidim_parscan_parvarcorrs(plot_type_flag,all_par_vals_lhs,stat_sol_nodes_lhs_parscan,...
                                    nodes,sel_nodes,[],[],[],plot_settings);

full_filename_with_path=fcn_save_fig(strjoin(plot_type_flag,'_'),save_folder,fig_file_type{1},'');                         

%% linear or lin-log regression of VARIABLES as fcn of PARAMETERS: VARIABLE=f(PARAMETER), the function plots R squared

plot_type_flag={'par_var','heatmap','r_sq'}; % {'par_var','heatmap'/'lineplot','r_sq'/'slope'}
sel_nodes=setdiff(3:numel(nodes),8); 
% plot_settings=[fontsize,maximum value for heatmap colors], if plot_settings(2)=NaN, then max color automatically selected
plot_settings=[14 NaN]; 
% if regression type is 'linlog', then the fit is y = a + b*log10(x)
regr_type={'log','linear'}; % linlog recommended if parameter values log-uniformly distributed in sampling
figure('name',strjoin(plot_type_flag))
[r_squared,slope_intercept]=fcn_multidim_parscan_parvarcorrs(plot_type_flag,all_par_vals_lhs,stat_sol_nodes_lhs_parscan,...
                                 nodes,sel_nodes,... % which nodes
                                 scan_params_sensit,scan_params_up_down_sensit, ... % parameters (CAREFUL that they are the same as in LHS!)
                                 regr_type{1},plot_settings);

fcn_save_fig(strjoin(plot_type_flag,'_'),save_folder,fig_file_type{1},'overwrite')

%% Quantify importance of parameters from LHS by a regression tree

% predictor importance values: look into arguments of <fitrtree> in MATLAB documentation to modulate regression tree
% for STATES or NODES?
scan_values=stat_sol_nodes_lhs_parscan; % stat_sol_states_lhs_parscan
sel_nodes=3:numel(nodes); % STATES or NODES to be analyzed
% names of selected transition rates and their predictor importance values
% [~,~,predictor_names] = fcn_get_trans_rates_tbl_inds(scan_params,scan_params_up_down,nodes);
% predictorImportance_vals=cell2mat(arrayfun(@(x) predictorImportance(...
%     fitrtree(all_par_vals_lhs,scan_values(:,x),'PredictorNames',predictor_names)), sel_nodes,'un',0)');

% CALCULATE and PLOT predictor importance
plot_type_flags={'line','bar'};
figure('name','regression_tree_pred_import')
[predictor_names,predictorImportance_vals] = fcn_multidim_parscan_predictorimport(scan_params_sensit,scan_params_up_down_sensit,...
                                                all_par_vals_lhs,scan_values,...
                                                nodes,sel_nodes,...
                                                plot_type_flags{2});

fcn_save_fig('regression_tree_pred_import',save_folder,fig_file_type{1},'')
                                            
%% Sobol total sensitivity metric                                            

% On Sobol total sensitivity index see: https://en.wikipedia.org/wiki/Variance-based_sensitivity_analysis
% This metric indicates how much of the total variance in a variable is due to variation in a given parameter
% We calculate here the usual numerical approximation of analytical equivalent from Monte Carlo sampling.
% From the LHS sampling above we take the matrices of parameter sets and variable values:
% [parameter sets, variable values]: [all_par_vals_lhs,stat_sol_nodes_lhs_parscan]

% select only parameters with an R-squared over some value
% [par_ind_table,~,~]=fcn_get_trans_rates_tbl_inds(scan_params_sensit,scan_params_up_down_sensit,nodes);
r_sq_thresh=0.2; par_ind_table_filtered=par_ind_table(sum(r_squared>r_sq_thresh)>0,:);
scan_params_filtered=unique(par_ind_table_filtered(:,1))'; % scan_params_unique=unique(par_ind_table(sum(r_squared>0.2)>0,1))'; 
scan_params_up_down_filtered=arrayfun(@(x) par_ind_table_filtered(par_ind_table_filtered(:,1)==x,2)', scan_params_filtered,'un',0);

% Sobol total sensitivity: calculated for one variable at a time
% selected nodes to display
sel_nodes=setdiff(1:numel(nodes),find(sum(cell2mat(arrayfun(@(x) strcmp(nodes,x), {'cc','KRAS','CDC25B'},'un',0)'))));
% sel_vars=[]; % if left empty, all nodes/states are analyzed
sample_size=300; % if left empty, the sample size is half of the original param scan <all_par_vals_lhs>
% PLOT SETTINGS: [fontsize_plot,fontsize_axes,fontsize_title, min_color(optional), max_color(opt), progress_calcul_every_x_% (opt)];
plot_settings=[14 14 22 NaN NaN 10];
var_types={'nodes','states'}; % analysis for states or nodes
% to calculate Sobol total sensitivity we need <sample_size*numel(scan_params_up_down)> evaluations of the model
figure('name','sobol sensitivity index')
sobol_sensit_index=fcn_multidim_parscan_sobol_sensit_index([],var_types{1},...
                            all_par_vals_lhs,stat_sol_nodes_lhs_parscan,stat_sol_states_lhs_parscan,...
                            sample_size,...
                            scan_params_filtered,scan_params_up_down_filtered,...% scan_params_sensit,scan_params_up_down_sensit
                            stg_table,x0,nodes,sel_nodes,plot_settings);

% if already calculated <sobol_sensit_index> and only want to plot results, provide <sobol_sensit_index> as FIRST argument 
fcn_multidim_parscan_sobol_sensit_index(sobol_sensit_index,var_types{1},[],[],[],[],...
                                scan_params_filtered,scan_params_up_down_filtered,[],[],nodes,sel_nodes,plot_settings);

% SAVE
fcn_save_fig('sobol_sensitivity_index',save_folder,fig_file_type{1},'y')

%% PARAMETER FITTING

% Relevant functions:
% transition_rates_table=fcn_trans_rates_table(nodes,'uniform',[],[],chosen_rates,chosen_rates_vals); % transition_rates_table=ones(size(transition_rates_table));
% tic; [A_sparse,~]=fcn_build_trans_matr(stg_table,transition_rates_table,''); toc; 
% tic; [stat_sol,term_verts_cell,cell_subgraphs]=split_calc_inverse(A_sparse,transition_rates_table,x0); toc
% [stationary_node_vals,init_node_vals]=fcn_calc_init_stat_nodevals(x0,stat_sol);

% define parameters to vary (predictor_names)
% sensitive parameters identified by 1-dimensional param scan: scan_params_sensit,scan_params_up_down_sensit
[~,~,predictor_names]=fcn_get_trans_rates_tbl_inds(scan_params_sensit,scan_params_up_down_sensit,nodes); 
% define data vector (generate some data OR load from elsewhere)
sel_param_vals=lognrnd(1,1,1,numel(predictor_names)); % abs(normrnd(1,0.5,1,numel(predictor_names)));
transition_rates_table=fcn_trans_rates_table(nodes,'uniform',[],[],predictor_names,sel_param_vals);
y_data=fcn_calc_init_stat_nodevals(x0,split_calc_inverse(fcn_build_trans_matr(stg_table,transition_rates_table,''),transition_rates_table,x0));

% create functions that calculate sum of squared deviations & values of
% variables (composed of different fcns) - RERUN THIS if you want to fit to new data or new non-fitted transition rates!!
[fcn_statsol_sum_sq_dev,fcn_statsol_values]=fcn_simul_anneal(y_data,x0,stg_table,nodes,predictor_names);

% FITTING by simulated annealing (look at arguments in anneal/anneal.m)
% initial guess for parameters
init_vals=rand(size(predictor_names)); init_error=fcn_statsol_sum_sq_dev(init_vals);
% simulated annealing with existing algorithm anneal/anneal.m, with modifications in script: 
% defined counter=0 before while loop, and inserted <T_loss(counter,:)=[T oldenergy];> at line 175, defined <T_loss> as 3rd output
% arguments for algorithm need to be provided as structure: 
% struct('Verbosity',2, 'StopVal', 0.01, 'StopTemp',1e-8) % stopping error value, stopping temperature parameter
fitting_arguments=struct('Verbosity',2, 'StopVal', 0.001);
tic; [optim_par_vals,best_error,T_loss]=anneal(fcn_statsol_sum_sq_dev,init_vals,fitting_arguments); toc 
% 20-40 mins for 15var KRAS model

% PLOT results
thres_ind=size(T_loss,1); % thres_ind=find(T_loss(:,2)<1e-2,1); 
vars_show=2; % 1=temperature, 2=error
figure('name','simul anneal')
plot(1:thres_ind, T_loss(1:thres_ind,vars_show),'LineWidth',4); xlim([0 thres_ind]); if init_error/best_error>30; set(gca,'yscale','log'); end
legend_strs = {'temperature', 'sum of squared error'}; legend(legend_strs{vars_show},'FontSize',22); grid on;
xlabel('number of iterations','FontSize',16); set(gca,'FontSize',16); title('Parameter fitting by simulated annealing','FontSize',22)
% SAVE FIGURE
fcn_save_fig(strcat('simulated_annealing_',num2str(numel(predictor_names)),'fittingpars'),save_folder,fig_file_type{1},'')

% mean absolute error: mean(abs(y_data - fcn_statsol_values(optim_par_vals)))
% distance of fitted params from true values: abs(optim_par_vals - sel_param_vals)./sel_param_vals