# exact-stoch-log-mod
Exact calculation of stationary states and parameter analysis of stochastic logical models

A MATLAB toolbox, version 2015b or later required.

[Wiki of the project](https://github.com/mbkoltai/exact-stoch-log-mod/wiki)

### To use the toolbox, the following freely available MATLAB toolboxes need to be downloaded and added to the path:
#### required for heatmaps:
[Customizable heatmaps](https://mathworks.com/matlabcentral/fileexchange/24253-customizable-heat-maps)  
addpath('heatmaps')

[Redblue colormap](https://mathworks.com/matlabcentral/fileexchange/25536-red-blue-colormap)  
addpath('redblue');

#### optional for figures with multiple subplots and to save figures
[tight subplots](https://mathworks.com/matlabcentral/fileexchange/27991-tight_subplot-nh-nw-gap-marg_h-marg_w) (for subplots with smaller gaps)  
addpath('tight_subplot') 

[export_fig](https://mathworks.com/matlabcentral/fileexchange/23629-export_fig) (export figures as EPS or PDF as they appear)  
addpath('altmany-export_fig-acfd348') 

[Simulated annealing](https://mathworks.com/matlabcentral/fileexchange/10548-general-simulated-annealing-algorithm) (parameter fitting by simulated annealing)  
addpath('anneal') 
