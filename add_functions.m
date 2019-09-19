
addpath('functions/') 
addpath('toolboxes/')
% using these publicly available MATLAB toolboxes:
% needed for PLOTS:
addpath('toolboxes/heatmaps') % https://mathworks.com/matlabcentral/fileexchange/24253-customizable-heat-maps
addpath('toolboxes/redblue'); % https://mathworks.com/matlabcentral/fileexchange/25536-red-blue-colormap (for redblue colormaps)

% optional for generating and saving plots
% for subplots with smaller gaps:
addpath('toolboxes/tight_subplot') % https://mathworks.com/matlabcentral/fileexchange/27991-tight_subplot-nh-nw-gap-marg_h-marg_w

% export figures as EPS or PDF as they appear on plots:
% Optional. https://mathworks.com/matlabcentral/fileexchange/23629-export_fig
export_fig_name=dir('toolboxes/altman*'); export_fig_name=export_fig_name.name; addpath(strcat('toolboxes/',export_fig_name)) 

% optional for PARAMETER FITTING by simulated annealing:
addpath('toolboxes/anneal') % https://mathworks.com/matlabcentral/fileexchange/10548-general-simulated-annealing-algorithm

% optional for maximal distinguishable colors on plots
addpath('toolboxes/distinguishable_colors/')

% optional defaults settings: title font weight normal, docked figures
set(0,'DefaultAxesTitleFontWeight','normal'); set(0,'DefaultFigureWindowStyle','docked');
