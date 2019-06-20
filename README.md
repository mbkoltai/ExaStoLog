# Exact-stoch-log-mod
Exact calculation of stationary states and parameter analysis of stochastic logical models

0. [Requirements](#markdown-header-0-toc)
1. [Model creation](#toc-model_creation)
2. [Calculation of stationary solution](#toc-calc-stat-sol)
3. [Visualizing the stationary solution](#toc-vis-stat-sol)
4. [Visualizing the state transition graph](#toc-vis-stg)
5. [Parameter sensitivity analysis: one-dimensional scans](#param-sens-1dim)
6. [Parameter sensitivity analysis: multi-dimensional scans](#param-sens-multidim)
7. [Parameter fitting by simulated annealing](#param-fitting)

### [Requirements](#markdown-header-0-toc)

#### - MATLAB version 2015b or later

#### - clone the repository and add the folder 'functions' to your path: addpath('functions')

#### - the following freely available MATLAB toolboxes also need to be downloaded and added to the path:

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


### [Model creation](#toc-model_creation)

### [Calculation of stationary solution](#toc-calc-stat-sol)

### [Visualizing the stationary solution](#toc-vis-stat-sol)

### [Visualizing the state transition graph](#toc-vis-stg)