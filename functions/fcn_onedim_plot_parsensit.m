function resp_coeff=fcn_onedim_plot_parsensit(plot_type_flag,var_type_flag,readout_type_flag,...
                                                scan_variable,nonzero_states_inds,parscan_matrix,...
                                                nodes,scan_params,scan_params_up_down,sensit_cutoff,param_settings)

% plot_type_flag: heatmap or lineplot
% var_type_flag: states or nodes
% readout_type_flag: variable value or response coefficient

[~,scan_par_inds,~]=fcn_get_trans_rates_tbl_inds(scan_params,scan_params_up_down,nodes);

fontsize_axes=param_settings(1);
fontsize_title=param_settings(2);

n_nodes=numel(nodes); truth_table_inputs=rem(floor([0:((2^n_nodes)-1)].'*pow2(0:-1:-n_nodes+1)),2);

trans_rates_names={strcat('u_',nodes);strcat('d_',nodes)}'; trans_rates_names=vertcat(trans_rates_names{:}); trans_rates_names=horzcat(trans_rates_names(:))';
nonzero_states=truth_table_inputs(nonzero_states_inds,:); 

% ONE STATE on one subplot a.a.f of all params
if strcmp(var_type_flag,'state') || strcmp(var_type_flag,'states')
    % scan_variable = stationary_state_vals_onedimscan;
    if size(scan_variable,3)~=size(nonzero_states,1)
        error('input variable (should be ''stationary_state_vals_onedimscan'') does not have correct dimension')
    end
elseif strcmp(var_type_flag,'node') || strcmp(var_type_flag,'nodes')
    % scan_variable = stationary_node_vals_onedimscan;
     if size(scan_variable,3)~=numel(nodes)
        error('input variable (should be ''stationary_node_vals_onedimscan'') does not have correct dimension')
    end
end

% calculate resp coeffs
num_diff_vals = diff(scan_variable,1,2); num_diff_parmatr = diff(parscan_matrix)';
% calc resp coeffs
resp_coeff=zeros(size(scan_variable)); resp_coeff = resp_coeff(:,2:end,:);
% resp_coeff_min_max = [floor(min(resp_coeff(:))*10)/10 ceil(max(resp_coeff(:))*10)/10];
for k=1:size(scan_variable,3)
    p_div_x = parscan_matrix'./scan_variable(:,:,k);
    resp_coeff(:,:,k) = (num_diff_vals(:,:,k)./num_diff_parmatr).*p_div_x(:,2:end);
end

% identify parameters that have an effect on stat vars
sensit_params_table=arrayfun(@(x) max(abs(resp_coeff(:,:,x)'))>sensit_cutoff, 1:size(resp_coeff,3),'un',0); sensit_params_table=vertcat(sensit_params_table{:});
sensit_pars=find(sum(sensit_params_table)>0);
sensit_vars=find(sum(sensit_params_table,2)>0)'; 

% take subspace of respcoeffs and statevars where there is an effect
resp_coeff_sensit_parts=resp_coeff(sensit_pars,:,sensit_vars);
scan_variable_sensit_parts=scan_variable(sum(sensit_params_table)>0,:,sum(sensit_params_table,2)>0);

% PLOTs showing the stationary value of variables
if strcmp(readout_type_flag,'var_value') || strcmp(readout_type_flag,'value') || strcmp(readout_type_flag,'values')

nrow=ceil(sqrt(size(scan_variable_sensit_parts,3))); ncol=ceil(sqrt(size(scan_variable_sensit_parts,3)));
% ONE NODE on one subplot a.a.f of all params

%%% LINEPLOT, VARIABLE VALUE 
if strcmp(plot_type_flag,'lineplot') || strcmp(plot_type_flag,'line')
    
nrow=ceil(sqrt(size(scan_variable,3))); ncol=ceil(sqrt(size(scan_variable,3)));

for k=1:size(scan_variable,3)
    subplot(nrow,ncol,k); 
    semilogx(parscan_matrix(:,sensit_pars), scan_variable(sensit_pars,:,k)', 'LineWidth',2); ylim([0 1]);
    
    if strcmp(var_type_flag,'node') || strcmp(var_type_flag,'nodes')
        title(strrep(nodes(k),'_','\_'), 'FontWeight','normal','FontSize',fontsize_title); ylim([0 1]); 
    else
        title(strcat('p([',num2str(nonzero_states(k,:)),'])'), 'FontWeight','normal','FontSize',fontsize_title); 
    end
    if rem(k,ncol)==1
        ylabel('stationary value','Fontsize',fontsize_axes)
    end
    
if k==size(scan_variable,3)
    legend(strrep(trans_rates_names(scan_par_inds(sensit_pars)),'_','\_'), 'Location', 'eastoutside');
end
end

%%%%%%%%%%%%%%%% 
%%% HEATMAP, VARIABLE VALUE
elseif strcmp(plot_type_flag,'heatmap') || strcmp(plot_type_flag,'hmap')
 
% resp_coeff_sensit_parts = reshape( resp_coeff_sensit_parts, size(resp_coeff_sensit_parts,2)*size(resp_coeff_sensit_parts,3), size(resp_coeff_sensit_parts,1) );
% resp_coeff_sensit_parts: [params, param values, variables]
for k=1:size(scan_variable_sensit_parts,3)
    if k>=size(scan_variable_sensit_parts,3)-ncol+1
        xlabs = trans_rates_names(scan_par_inds(sum(sensit_params_table)>0));
    else
        xlabs =[];
    end
subplot(nrow,ncol,k);    
var_to_plot = scan_variable_sensit_parts(:,:,k)'; % resp_coeff_sensit_parts
heatmap( var_to_plot,xlabs,[],[],'TickAngle',90,'Colormap','redblue',...
    'MinColorValue',0,'MaxColorValue',1,'GridLines','none','FontSize',11,'ShowAllTicks',true);    
    % max(abs(resp_coeff_min_max))
    
    if strcmp(var_type_flag,'node') || strcmp(var_type_flag,'nodes')
        title(strrep(nodes(sensit_vars(k)),'_','\_'), 'FontWeight','normal','FontSize',fontsize_title)
    else
        title(strcat('p([',num2str(nonzero_states(k,:)),'])'), 'FontWeight','normal','FontSize',fontsize_title);
    end
    
    if k==size(scan_variable_sensit_parts,3)
        colorbar('EastOutside');
    end
    
    if rem(k,ncol)==1
        ylabel('stationary value', 'FontSize', fontsize_axes)
    end

    
end % end of for loop

    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PLOT RESPONSE COEFFICIENTS (sensitivity metric)

elseif strcmp(readout_type_flag,'sensitivity') || strcmp(readout_type_flag,'respcoeff') || strcmp(readout_type_flag,'resp_coeff')

nrow=ceil(sqrt(size(resp_coeff_sensit_parts,3))); ncol=ceil(sqrt(size(resp_coeff_sensit_parts,3)));

%%%%%%%%%%%%    
% RESP COEFF, LINEPLOT
if strcmp(plot_type_flag,'lineplot') || strcmp(plot_type_flag,'line')
% LINEPLOT w resp coeffs
% sensit_cutoff=0.04;
for k=1:size(resp_coeff_sensit_parts,3)
%     p_div_x = parscan_matrix'./stationary_node_vals_onedimscan(:,:,k);
%     resp_coeff(:,:,k) = (num_diff_vals(:,:,k)./num_diff_parmatr).*p_div_x(:,2:end);
% PLOT
subplot(nrow,ncol,k)
resp_coeff_var=resp_coeff_sensit_parts(:,:,k)'; % sensit_pars=max(abs(resp_coeff_var))>sensit_cutoff;
% lineplot
semilogx(parscan_matrix(2:end,sensit_pars), resp_coeff_var, 'LineWidth',2); 
ylim([floor(min(resp_coeff_sensit_parts(:))*10)/10 ceil(max(resp_coeff_sensit_parts(:))*10)/10])
if k==size(resp_coeff_sensit_parts,3) 
    legend(strrep(trans_rates_names(scan_par_inds(sensit_pars)),'_','\_'), 'Location', 'EastOutside');
end
if rem(k,ncol)==1
    ylabel('response coefficient', 'FontSize', fontsize_axes)
end

if strcmp(var_type_flag,'node') || strcmp(var_type_flag,'nodes')
        title(strrep(nodes(sensit_vars(k)),'_','\_'), 'FontWeight','normal','FontSize',fontsize_title)
else
        title(strcat('p([',num2str(nonzero_states(k,:)),'])'), 'FontWeight','normal','FontSize',fontsize_title); 
end

end

%%%%%%%%%%%%    
% RESP COEFF, HEATMAP
else
% heatmap
% parameters that have abs(resp_coeff)>cutoff
sensit_params_table=arrayfun(@(x) max(abs(resp_coeff(:,:,x)'))>sensit_cutoff, 1:size(resp_coeff,3),'un',0); sensit_params_table=vertcat(sensit_params_table{:});
resp_coeff_sensit_parts=resp_coeff(sum(sensit_params_table)>0,:,sum(sensit_params_table,2)>0);

sensit_vars=find(sum(sensit_params_table,2)>0)';
% resp_coeff_sensit_parts = reshape( resp_coeff_sensit_parts, size(resp_coeff_sensit_parts,2)*size(resp_coeff_sensit_parts,3), size(resp_coeff_sensit_parts,1) );
% resp_coeff_sensit_parts: [params, param values, variables]
color_min_max=[floor(min(resp_coeff_sensit_parts(:))*10)/10 ceil(max(resp_coeff_sensit_parts(:))*10)/10];
for k=1:size(resp_coeff_sensit_parts,3)
    if k>=size(resp_coeff_sensit_parts,3)-ncol+1
        xlabs = trans_rates_names(scan_par_inds(sum(sensit_params_table)>0));
    else
        xlabs =[];
    end
subplot(nrow,ncol,k);    
var_to_plot = resp_coeff_sensit_parts(:,:,k)'; % resp_coeff_sensit_parts

heatmap( var_to_plot,xlabs,[],[],'TickAngle',90,'Colormap','redblue',...
    'MinColorValue',color_min_max(1),'MaxColorValue',color_min_max(2),'GridLines','none','FontSize',11,'ShowAllTicks',true);    
    % max(abs(resp_coeff_min_max))
    
    if strcmp(var_type_flag,'node') || strcmp(var_type_flag,'nodes')
        title(strrep(nodes(sensit_vars(k)),'_','\_'), 'FontWeight','normal','FontSize',fontsize_title)
    else
        title(strcat('p([',num2str(nonzero_states(k,:)),'])'), 'FontWeight','normal','FontSize',fontsize_title); 
    end

    
    if rem(k,ncol)==1
        ylabel('response coeff.','FontSize',fontsize_axes)
    end
    if k==size(resp_coeff_sensit_parts,3)
        colorbar('EastOutside');
    end
    
end

end % plot type

else % readout type (value or sensitivity)
   
    error('readout_type_flag should be ''values/var_value'' or ''respcoeff/sensitivity'' ')
    
end