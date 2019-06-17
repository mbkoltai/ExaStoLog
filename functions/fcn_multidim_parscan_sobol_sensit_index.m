function tau_i=fcn_multidim_parscan_sobol_sensit_index(sobol_sensit_index,all_par_vals_lhs,stat_sol_lhs_parscan,sample_size,...
                                        scan_params,scan_params_up_down,stg_table,x0,nodes,sel_nodes,plot_settings)

par_ind_table=[repelem(scan_params, cellfun(@(x) numel(x),scan_params_up_down))', horzcat(scan_params_up_down{:})'];
prefix={'u_','d_'}; predictor_names=arrayfun(@(x) strcat(prefix(par_ind_table(x,2)), nodes(par_ind_table(x,1))), 1:size(par_ind_table,1),'un',0);
predictor_names=vertcat(predictor_names{:})';
                                    
if isempty(sobol_sensit_index)

if ~isempty(sample_size)
M=sample_size;
else
M=size(all_par_vals_lhs,1)/2; 
end

A=all_par_vals_lhs(1:M,:); B=all_par_vals_lhs(M+1:2*M,:);
% f_B_i_store=zeros(M,numel(nodes),size(all_par_vals_lhs,2)); 
tau_i=zeros(size(all_par_vals_lhs,2),numel(sel_nodes));
transition_rates_table=ones(2,numel(nodes));

for k = 1:size(all_par_vals_lhs,2)
    %%%%
    B_i = A; B_i(:,k) = B(:,k);
    % rerun calculations
    f_B_i=fcn_calc_paramsample_table(B_i,scan_params,scan_params_up_down,transition_rates_table,stg_table,x0,nodes);
    % f_B_i_store(:,:,k) = f_B_i;
    % sensit index
for varcount=1:numel(sel_nodes)
    diff_fA_fBi=( stat_sol_lhs_parscan(1:M,sel_nodes(varcount))-f_B_i(:,sel_nodes(varcount)) );
	tau_i(k,varcount)=(diff_fA_fBi'*diff_fA_fBi)/( 2*M*var(stat_sol_lhs_parscan(1:M,sel_nodes(varcount))) );
end
disp(k)
end

else
    tau_i=sobol_sensit_index;
end

num_size_plot=plot_settings(1); min_col_val=0; maxval_color=1; % 1.05*max(tau_i(:))
heatmap(tau_i,nodes(sel_nodes),predictor_names,'%0.2f','TickAngle',90,'Colormap','redblue',...
    'MinColorValue',min_col_val,'MaxColorValue',maxval_color,'GridLines','-','FontSize',num_size_plot,'ShowAllTicks',true,'colorbar',true)
set(gca,'FontSize',plot_settings(2)); title('Sobol total sensitivity index', 'Fontweight','normal', 'FontSize',plot_settings(3));
