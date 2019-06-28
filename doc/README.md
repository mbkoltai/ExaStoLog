# Exact-stoch-log-mod

Exact calculation of stationary states + parameter analysis & fitting of stochastic logical models

## Table of contents

1. [Requirements](#1-requirements)
1. [Model creation](#2-model-creation)
1. [Calculation of stationary solution](#3-calculation-of-stationary-solution)
1. [Visualizing the stationary solution](#4-visualizing-the-stationary-solution)
1. [Visualizing the state transition graph](#5-visualizing-the-state-transition-graph)
1. [One-dimensional parameter sensitivity analysis](#6-one-dimensional-parameter-sensitivity-analysis)
1. [Multi-dimensional parameter sensitivity analysis](#7-multi-dimensional-parameter-sensitivity-analysis)
1. [Parameter fitting by simulated annealing](#8-parameter-fitting-by-simulated-annealing)

### 1. Prerequisites

##### - MATLAB version 2015b or later

##### - clone the [repository](https://github.com/mbkoltai/exact-stoch-log-mod) and add the folder 'functions' to your path: addpath('functions')

##### - the following freely available MATLAB toolboxes also need to be downloaded and added to the path:

- [Customizable heatmaps](https://mathworks.com/matlabcentral/fileexchange/24253-customizable-heat-maps)  
addpath('heatmaps')

- [Redblue colormap](https://mathworks.com/matlabcentral/fileexchange/25536-red-blue-colormap)  
addpath('redblue');

Optional (for figures with multiple subplots, to save figures and for parameter fitting):  
- [tight subplots](https://mathworks.com/matlabcentral/fileexchange/27991-tight_subplot-nh-nw-gap-marg_h-marg_w) (for subplots with smaller gaps)  
addpath('tight_subplot')

- [export_fig](https://mathworks.com/matlabcentral/fileexchange/23629-export_fig) (export figures as EPS or PDF as they appear)  
addpath('altmany-export_fig-acfd348')

- [Simulated annealing](https://mathworks.com/matlabcentral/fileexchange/10548-general-simulated-annealing-algorithm) (parameter fitting by simulated annealing)  
addpath('anneal')  
To plot the convergence of the fitting process modify the script by 1) defining \<T_loss\> as 3rd output of the function, 2) inserting <counter=0> before the while loop
and 3) inserting <T_loss(counter,:)=[T oldenergy];> at line 175 within the while loop.


### 2. Model creation

Models can be defined by entering the list of nodes and their corresponding rules as a cell of strings, using MATLAB logical notation ('&', '|', '~', '(', ')'), for instance the following is a 10-node model of mitotic entry:  
```MATLAB
nodes = {'cc','kras', 'dna_dam', 'chek1', 'mk2', 'atm_atr', 'hr','cdc25b', 'g2m_trans', 'cell_death'};

rules={'cc',...  
'kras',...  
'(dna_dam | kras) & ~hr',...  
'atm_atr',...  
'atm_atr & kras',...  
'dna_dam',...  
'(atm_atr  | hr) & ~cell_death',...  
'(cc|kras) & (~chek1 & ~mk2) & ~cell_death',...  
'g2m_trans | cdc25b',...  
'cell_death | (dna_dam & g2m_trans)'};  
```


Alternatively, models can be read in in boolnet (.bnet) format by defining the model's path and reading in the boolnet file by the custom function *fcn_bnet_readin*:  
```MATLAB
% name of the model
model_name='mammalian_cc'; % kras15vars

% model read in from an existing BOOLNET file
[nodes,rules]=fcn_bnet_readin('model_files/traynard2016_mammalian_cellcycle.bnet'); krasmodel10vars.bnet
```

Once we have the list of nodes and their logical rules, we can check if all variables referred to by rules are found in the list of nodes:
```MATLAB
fcn_nodes_rules_cmp(nodes,rules)
```

If the BOOLNET file is not consistent in terms of its nodes and rules, this function will display an error message, otherwise prints *Model seems correct: all elements in rules found in nodes list*.

To create the logical model we need to generate a function file:
```MATLAB
truth_table_filename='fcn_truthtable.m';
fcn_write_logicrules(nodes,rules,truth_table_filename)
```
**State transition graph**

From this function file we generate the state transition graph (STG) of the logical model, that is still independent of values of transition rates, it only informs us amongst which states we have possible transitions in the system:
```MATLAB
[stg_table,~,~]=fcn_build_stg_table(truth_table_filename,nodes,'','');
```
This step does not need to be repeated as the STG does not change with the parameter values, eg. for parameter sensitivity analysis.

Below is the plot for the STG of the aforementioned 10-node model (the function for this plot is described in [Visualizing the state transition graph](#5-visualizing-the-state-transition-graph)):
![STG_10nodes_full_subgraph](../sample_plots/kras10vars/STG_10nodes_full_subgraph.png)




Since this model has 2 inputs whose state cannot change, its global STG is composed of 4 disconnected subgraphs, one of which is shown on the right.

### 3. Calculation of stationary solution

**Defining transition rates**

To calculate the stationary states of a model we need to assign values to the *2xn* (n=number of nodes) transition rates of the model.

We can select a subset of the transition rates to have a specific value by their names, which is always made up of *d_* or *u_* and the name of the respective node, and also define the vector of values we want them to have:
```MATLAB
chosen_rates={'u_cdc25b','d_dna_dam'}; chosen_rates_vals=[0.25, 0.15];
```

Otherwise, if you want all transition rates to have the same value or to be sampled from the same random distribution, leave these variables empty (or don't define them):
```MATLAB
chosen_rates=[]; chosen_rates_vals=[];
```
Call the function to generate the table of transition rates, selecting if you want to have the rates to have a uniform value or to be sampled from a normal distribution, in this case specify the mean and standard distribution:
```MATLAB
% ARGUMENTS
% <uniform> assigns a value of 1 to all params. other option: <random>
distr_type={'uniform','random'};
% if 'random' is chosen, the mean and standard dev of a normal distrib has to be defined
meanval=[]; sd_val=[];
transition_rates_table=fcn_trans_rates_table(nodes,distr_type{1},meanval,sd_val,chosen_rates,chosen_rates_vals);
```

**Creating transition matrix**

Now we can build the transition matrix of the model with the specified transition rates:
```MATLAB
[A_sparse,~]=fcn_build_trans_matr(stg_table,transition_rates_table,'');
```

For the subsequent calculations only the transition matrix *A* is needed (it is converted to the kinetic matrix within functions), if you want to have the kinetic matrix *K* as a variable (*dp(t)/dt=Kp(t)*, as opposed to *p(t+1)=p(t)A*), then run the function as:
```MATLAB
[A_sparse,K_sparse]=fcn_build_trans_matr(stg_table,transition_rates_table,'kinetic');
```

**Defining initial conditions**

We can define an initial condition by specifying a given state we want to have a larger than random (there are 2^n states in total) probability. To do this select the nodes you want to have a value of 1 at the initial state by their names and the probability value of this state:
```MATLAB
initial_on_nodes = {'CycD','Rb_b1','Rb_b2','Cdh1','p27_b1','p27_b2','Skp2'};
% what is the probability of this state, (eg. dom_prob=0.8, ie. 80% probability)
dom_prob=0.8;
```

Then with the function *fcn_define_initial_states* we assing the remaining *dom_prob* probability either among all other states (*distrib_types='broad'*) or all states where the selected nodes have a value of 1 (*distrib_types='restrict'*):
```MATLAB
distrib_types={'restrict','broad'}; plot_flag=[]; % if plot_flag non-empty, we get a bar plot of initial values
x0=fcn_define_initial_states(initial_on_nodes,dom_prob,nodes,distrib_types{1},plot_flag);
```
If \<plot_flag> is non-empty we also get a bar plot of the initial states.

You can also define a completely random probability distribution of initial states:
```MATLAB
n_nodes=numel(nodes); truth_table_inputs=rem(floor([0:((2^n_nodes)-1)].'*pow2(0:-1:-n_nodes+1)),2);
x0=zeros(1,2^n_nodes)'; x0=rand(1,size(truth_table_inputs,1))'; x0=x0/sum(x0);
```

**Calculation of stationary solution**

With the transition matrix, table of transition rates and the initial condition defined we can now calculate the stationary solution of the model (it is informative to time this calculation for larger calculations later on when it will be repeated):

```MATLAB
tic; [stat_sol,term_verts_cell,cell_subgraphs]=split_calc_inverse(A_sparse,transition_rates_table,x0); toc
```

The outputs of the calculation are:  
**stat_sol**: stationary solution for all the states  
**term_verts_cell**: index of nonzero states. If the STG is disconnected, then the nonzero states corresp to these disconn subgraphs are in separate cells  
**cell_subgraphs**: indices of states belonging to disconnected subgraphs (if any)

The nonzero states can be quickly queried by:
```MATLAB
stat_sol(stat_sol>0)' % probability values of nonzero states
truth_table_inputs(stat_sol>0,:) % logical states that are nonzero
```

To have the stationary solution (and also the initial states) in terms of the probabilities of the model's _nodes_ having a value of 1, call the function:
```MATLAB
[stationary_node_vals,init_node_vals]=fcn_calc_init_stat_nodevals(x0,stat_sol);
```

### 4. Visualizing the stationary solution

We can now visualize the stationary solution, along with the transition (or kinetic) matrix of the model. Visualization of the full transition matrix can take very long, it is not recommended for models larger than 12 nodes.
The following arguments need to be defined for the visualization:
```MATLAB
% ARGUMENTS
% matrix_input: [], K_sparse or A_sparse (kinetic or transition matrix)
% min_max_col: minimum and max color for heatmap
% fontsize: [font size on the heatmap, title font size for stationary solutions]
% barwidth_states_val: width of the bars for bar plot of stationary solutions of states
% sel_nodes: nodes to show. If left empty, all nodes are shown
% fontsize_hm: fonsize on heatmap,
% fontsize_stat_sol: fontsize on barplots for stationary solution
% nonzero_flag: minimal value for probability to display - if this is non-empty, only plot nonzero states, useful for visibility if there are many states
```

Call the function with matrix visualization:
```MATLAB
sel_nodes=[]; min_max_col=[0 1]; barwidth_states_val=0.8; fontsize_hm=10; fontsize_stat_sol=20;
fontsize=[fontsize_hm fontsize_stat_sol]; matrix_input=A_sparse; plot_settings = [fontsize barwidth_states_val min_max_col]; nonzero_flag=0.01;

fcn_plot_A_K_stat_sol(matrix_input, nodes, sel_nodes, stat_sol, x0, plot_settings ,nonzero_flag)
```

![kras10vars_K_statsol](../sample_plots/kras10vars/K_statsol_10vars.png)

Or, without the transition matrix, with only the stationary solutions:

```MATLAB
sel_nodes=[]; nonzero_flag=0.01; barwidth_states_val=0.8; % for 3 nonzero states ~0.8 is a good value
fontsize=[9 20]; % [fontsize_y_axis_states,fontsize_x_axes_and_titles]
plot_settings=[fontsize barwidth_states_val]; matrix_input=[];
fcn_plot_A_K_stat_sol(matrix_input, nodes, sel_nodes, stat_sol, x0, plot_settings,nonzero_flag)
```

![kras15vars_statsol](../sample_plots/kras15vars/statsol_15vars.png)

Save the plot by running (**export_fig** toolbox needed):
```MATLAB
if exist(strcat(save_folder,model_name),'dir')==0; mkdir(strcat(save_folder,model_name)); end
fig_file_type={'.png','.eps'}; if ~isempty(matrix_input); matrix_input_str='_with_matrix'; else matrix_input_str=''; end
export_fig(strcat(save_folder,model_name,'/','single_solution_states_nodes_stat_sol',matrix_input_str,fig_file_type{2}),'-transparent','-nocrop')
```

**Visualize binary heatmap of nonzero stationary states**

To visualize what are the states that are the fixed points or cyclic attractor of your model on a heatmap, provide the following arguments and call the plotting function:
```MATLAB
% term_verts_cell: which subgraph to plot if there are disconnected ~
% num_size_plot: font size of 0/1s on the heatmap
% hor_gap: horizontal gap between terminal SCCs, bottom_marg: bottom margin, left_marg: left margin
numsize_plot=12; fontsize=12; hor_gap=0.02; bottom_marg=0.1; left_marg=0.04;
param_settings=[numsize_plot fontsize hor_gap bottom_marg left_marg];
% index of nonempty subgraph, check by <term_verts_cell>
nonempty_subgraph=2;
% want to use tight subplot? | order states by probability?
tight_subplot_flag='yes'; ranking_flag='yes';
% nodes to show. if none selected, then all nodes shown
sel_nodes=2:numel(nodes)-1;
% probability threshold for states to show (if left empty, all states shown)
prob_thresh=0.01;  % []; % 0.05;
% PLOT
statsol_binary_heatmap=fcn_plot_statsol_bin_hmap(stat_sol,prob_thresh,term_verts_cell{nonempty_subgraph},...
                            nodes,sel_nodes,param_settings,tight_subplot_flag,ranking_flag);
```

![binary_statsol_heatmap_15vars](../sample_plots/kras15vars/binary_statsol_heatmap_15vars.png)

The plot

### 5. Visualizing the state transition graph

### 6. One-dimensional parameter sensitivity analysis

### 7. Multi-dimensional parameter sensitivity analysis

![LHS_multidimparscankras10vars](readmeplots/LHS_multidimparscankras10vars.png)

### 8. Parameter fitting by simulated annealing
