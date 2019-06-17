function [scan_pars_sensit,scan_params_sensit_up_down]=fcn_onedim_sensit_pars_resp_coeff(resp_coeff,sensit_cutoff,...
    scan_params,scan_params_up_down,nodes)

sensit_params_table=arrayfun(@(x) max(abs(resp_coeff(:,:,x)'))>sensit_cutoff, 1:size(resp_coeff,3),'un',0); sensit_params_table=vertcat(sensit_params_table{:});
sensit_pars=find(sum(sensit_params_table)>0);
[scan_par_table,~]= fcn_get_trans_rates_tbl_inds(scan_params,scan_params_up_down,nodes);
% sensitive parameters
scan_pars_sensit=unique(scan_par_table(sensit_pars,1))';
for k=1:numel(scan_pars_sensit)
    scan_params_sensit_up_down{k}=scan_par_table(sensit_pars(scan_par_table(sensit_pars,1)==scan_pars_sensit(k)),2)';
    % sensit_par_lengths=cell2mat(arrayfun(@(x) numel(scan_params_up_down_sensit{x}), 1:numel(scan_params_up_down_sensit),'un',0));
end
