# ExaStoLog

A MATLAB toolbox for the exact calculation of stationary states + parameter sensitivity analysis & fitting of stochastic logical models.  
Author: Mihály Koltai, [Computational Systems Biology of Cancer at Institut Curie](https://github.com/sysbio-curie)

Read the tutorial [here](https://github.com/mbkoltai/exact-stoch-log-mod/tree/master/doc) (**under construction**).

### Requirements

#### - MATLAB version 2015b or later

#### - clone the repository and add the folder 'functions' to your path: addpath('functions')

#### - the following freely available MATLAB toolboxes also need to be downloaded and added to the path:

- [Customizable heatmaps](https://mathworks.com/matlabcentral/fileexchange/24253-customizable-heat-maps)  
addpath('heatmaps')

- [Redblue colormap](https://mathworks.com/matlabcentral/fileexchange/25536-red-blue-colormap)  
addpath('redblue');

Optional (for figures with multiple subplots, to save high-quality figures and for parameter fitting):  
- [tight subplots](https://mathworks.com/matlabcentral/fileexchange/27991-tight_subplot-nh-nw-gap-marg_h-marg_w) (for subplots with smaller gaps)  
addpath('tight_subplot') 

- [export_fig](https://mathworks.com/matlabcentral/fileexchange/23629-export_fig) (export figures as EPS or PDF as they appear)  
addpath('altmany-export_fig-acfd348') 

- [Simulated annealing](https://mathworks.com/matlabcentral/fileexchange/10548-general-simulated-annealing-algorithm) (parameter fitting by simulated annealing)  
addpath('anneal') 
