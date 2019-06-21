# Exact-stoch-log-mod

Exact calculation of stationary states + parameter analysis & fitting of stochastic logical models

## Table of contents

1. [Requirements](#1-requirements)
1. [Model creation](#2-model-creation)
1. [Calculation of stationary solution](#3-calc-stat-sol)
1. [Visualizing the stationary solution](#4-vis-stat-sol)
1. [Visualizing the state transition graph](#5-vis-stg)
1. [Parameter sensitivity analysis: one-dimensional scans](#6-sens-1dim)
1. [Parameter sensitivity analysis: multi-dimensional scans](#7-sens-multidim)
1. [Parameter fitting by simulated annealing](#8-param-fitting)

### 1. Requirements

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


### 2. Model creation

### 3. Calculation of stationary solution

### 4. Visualizing the stationary solution

### 5. Visualizing the state transition graph

### 6. Parameter sensitivity analysis: one-dimensional scans

### 7. Parameter sensitivity analysis: multi-dimensional scans

### 8. Parameter fitting by simulated annealing
