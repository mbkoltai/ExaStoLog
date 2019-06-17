function [state_transitions_inds,K_sparse,A_sparse_fast]=fcn_build_stg_table(truth_table_filename,nodes,transition_rates_table,num_matrix_flag)

n_nodes=numel(nodes); all_binary_states=rem(floor([0:((2^n_nodes)-1)].'*pow2(0:-1:-n_nodes+1)),2); dim_matr=2^n_nodes;

% delete(fcn_name)
if isa(strrep(truth_table_filename,'.m',''),'function_handle')
    disp('handle exists!!')
    clear(fcn_name) %  clear(func2str(yourfunctionhandle))
end
fcn_name_str=strrep(truth_table_filename,'.m',''); fcn_name=str2func(fcn_name_str); 
update_table=fcn_name(nodes); 

% if size(update_table,2)~=n_nodes
%     clear(fcn_name);
%     update_table=fcn_name(nodes);
% end

trans_source_states=cell(n_nodes,1); trans_target_states=trans_source_states; 
up_down_inds=trans_source_states; node_inds=up_down_inds; update_asynchr_table=zeros(2^n_nodes,n_nodes);

for k=1:n_nodes
    % disp(strcat('building ', num2str(k),'th variable transitions'));
    update_asynchr_table(:,setdiff(1:n_nodes,k))=all_binary_states(:,setdiff(1:n_nodes,k)); 
    update_asynchr_table(:,k)=update_table(:,k);
    trans_source_states{k}=find(all_binary_states(:,k)~=update_asynchr_table(:,k)); 
    [~,I_outputs,I_inputs]=intersect(all_binary_states,update_asynchr_table(trans_source_states{k},:),'rows'); 
    % update_asynchr_table(row_inds(trans_source_states{k}),:)
    n_trans = numel(I_inputs); [~,b]=intersect(I_inputs,(1:n_trans)');
    trans_target_states{k} = I_outputs(b);
    up_down_inds{k} = sum(update_asynchr_table(trans_source_states{k},:)-all_binary_states(trans_source_states{k},:),2);
    node_inds{k} = repelem(k,n_trans,1);
end

% convert -1 to 2 (down rates)
up_down_inds_mat=cell2mat(up_down_inds); trans_source_states_mat=cell2mat(trans_source_states); trans_target_states_mat=cell2mat(trans_target_states);
up_down_inds_arr = sign(up_down_inds_mat).^2 + (up_down_inds_mat<0);

% [source_state target_state node_ind up_down_ind]
state_transitions_inds=[trans_source_states_mat trans_target_states_mat cell2mat(node_inds) up_down_inds_arr];
    
if ~isempty(num_matrix_flag)
rate_inds=sub2ind(size(transition_rates_table),up_down_inds_arr,cell2mat(node_inds));

A_sparse_fast=sparse(trans_source_states_mat,trans_target_states_mat, ...
    transition_rates_table(rate_inds)/sum(transition_rates_table(:)),dim_matr,dim_matr);
diag_vals=1-sum(A_sparse_fast,2); dim_matr_arr=(1:dim_matr)';
rows = [trans_source_states_mat;dim_matr_arr]; cols=[trans_target_states_mat;dim_matr_arr]; 
vals=[transition_rates_table(rate_inds)/sum(transition_rates_table(:));diag_vals];
A_sparse_fast = sparse(rows,cols,vals,dim_matr,dim_matr);
K_sparse=(transpose(A_sparse_fast) - speye(size(A_sparse_fast)) )*sum(transition_rates_table(:));

else
    A_sparse_fast=[]; K_sparse=[];
end

