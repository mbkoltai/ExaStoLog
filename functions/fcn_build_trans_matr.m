function [K_sparse,A_sparse_fast]=fcn_build_trans_matr(state_transitions_inds,transition_rates_table,nodes,kin_matr_flag)

dim_matr=2^numel(nodes);

% state_transitions_inds=[trans_source_states_mat, trans_target_states_mat, cell2mat(node_inds), up_down_inds_arr];
    
trans_source_states_mat=state_transitions_inds(:,1); 
trans_target_states_mat=state_transitions_inds(:,2); 
up_down_inds_arr=state_transitions_inds(:,4); 

rate_inds=sub2ind(size(transition_rates_table),up_down_inds_arr, state_transitions_inds(:,3));

A_sparse_fast=sparse(trans_source_states_mat,trans_target_states_mat, ...
    transition_rates_table(rate_inds)/sum(transition_rates_table(:)),dim_matr,dim_matr);
diag_vals=1-sum(A_sparse_fast,2); dim_matr_arr=(1:dim_matr)';
rows = [trans_source_states_mat;dim_matr_arr]; cols=[trans_target_states_mat;dim_matr_arr]; 
vals=[transition_rates_table(rate_inds)/sum(transition_rates_table(:));diag_vals];
A_sparse_fast = sparse(rows,cols,vals,dim_matr,dim_matr);
if ~isempty(kin_matr_flag)
    K_sparse=(transpose(A_sparse_fast) - speye(size(A_sparse_fast)) )*sum(transition_rates_table(:));
else
    K_sparse=[];
end
