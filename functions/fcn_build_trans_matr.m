function [A_sparse_fast,K_sparse]=fcn_build_trans_matr(stg_table,transition_rates_table,kin_matr_flag)

dim_matr=2^size(transition_rates_table,2);

% state_transitions_inds=[trans_source_states_mat, trans_target_states_mat, cell2mat(node_inds), up_down_inds_arr];
    
% trans_source_states_mat=stg_table(:,1); trans_target_states_mat=stg_table(:,2); up_down_inds_arr=stg_table(:,4); 
rate_inds=(stg_table(:,3)-1)*2+stg_table(:,4);
% sub2ind(size(transition_rates_table),up_down_inds_arr, stg_table(:,3));

B=sparse(stg_table(:,1), stg_table(:,2), transition_rates_table(rate_inds)/sum(transition_rates_table(:)),dim_matr,dim_matr);
% diag_vals=1-sum(A_sparse_fast,2); % dim_matr_arr=(1:dim_matr)';
% rows = [stg_table(:,1); (1:dim_matr)']; 
% cols=[stg_table(:,2); (1:dim_matr)']; 
% vals=[transition_rates_table(rate_inds)/sum(transition_rates_table(:));diag_vals];
% A_sparse_fast = sparse(rows,cols,vals,dim_matr,dim_matr);
A_sparse_fast = B + (speye(size(B)) - diag(sum(B,2)));

if ~isempty(kin_matr_flag)
    K_sparse=(transpose(A_sparse_fast) - speye(size(A_sparse_fast)) )*sum(transition_rates_table(:));
else
    K_sparse=[];
end
