function [stationary_node_vals,init_node_vals]=fcn_calc_init_stat_nodevals(x0,stat_sol)

init_node_vals=[];
n=log2(size(x0,1)); truth_table_inputs=rem(floor([0:((2^n)-1)].'*pow2(0:-1:-n+1)),2); 
if (numel(stat_sol)==numel(x0)) && (size(stat_sol,1)==size(x0,1)) && 2^n==numel(stat_sol)
stationary_node_vals=stat_sol'*truth_table_inputs; 
if ~isempty(x0)
    init_node_vals=x0'*truth_table_inputs;
end

else
   disp('initial condition or stationary solution don''t match with each other and/or with number of nodes')
end