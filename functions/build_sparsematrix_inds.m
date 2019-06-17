function [state_transition_table,K_sparse_table,K_sparse,A_sparse]=build_sparsematrix_inds(update_matrix,nodes,transition_rates_table)

n=numel(nodes);
truth_table_inputs = rem(floor([0:((2^n)-1)].'* pow2(0:-1:-n+1)),2);

state_no=size(truth_table_inputs,1); l_nodes=length(nodes);

% number of states and nodes         
if size(truth_table_inputs,2)==l_nodes

% maximal dimension of transition matrix (=rows of transition table) is states^2
if (state_no^2)<2e6
state_transition_table=repelem(nan,state_no^2,4);
end

for k = 1:length(nodes)
    % disp(nodes(k));
    % the update matrix has all nodes updated
    % here we update per node, so output_table has same columns as truth_table_inputs, except in column 'i'
    output_table = zeros(state_no,l_nodes); output_table(:,k)=update_matrix(:,k); output_table(:,setdiff(1:l_nodes,k))=truth_table_inputs(:,setdiff(1:l_nodes,k));
    % take only those rows where there is a change, input state=/=output state
    % these are also the row indices of 'truth_table_inputs' where we have output states different from source states
    nonzero_transitions_ind = find(sum(truth_table_inputs==output_table,2)<l_nodes);
    poss_trans_table_outputs=output_table(nonzero_transitions_ind ,:);
    % poss_trans_table_inputs=truth_table_inputs(sum(truth_table_inputs==output_table,2)<l_nodes,:);
    
    % find index of output states, IA will have indices from input states
    [~, I_outputs, I_inputs] = intersect(truth_table_inputs,poss_trans_table_outputs,'rows');
    % [~, I_inputs] = intersect(truth_table_inputs,poss_trans_table_inputs,'rows');
    % indices of rows from output table
    truth_table_rows=find(sum(truth_table_inputs==output_table,2)<l_nodes);
    
    % truth_table_inputs(truth_table_rows,:) poss_trans_table 
    nrow=size(truth_table_rows,1);
    
    if k==1
        nrow_cumul=nrow;
    else
        nrow_cumul=nrow_cumul+nrow;
    end
    
	if nrow>0
        % first value is the node index, second is whether up (1) or down (2) rate
        param_node_updown_ind = [repelem(k,numel(I_inputs),1) truth_table_inputs(nonzero_transitions_ind(I_inputs),k)+1];
        % sub2ind(param_table_size,param_node_updown_ind(2,:),param_node_updown_ind(1,:))
      state_transition_table((nrow_cumul-nrow+1):nrow_cumul,:) = [nonzero_transitions_ind(I_inputs) I_outputs param_node_updown_ind(:,1) param_node_updown_ind(:,2)];
	end
    % [trans_matrix_indices(:,1:2) sub2ind(param_table_size,trans_matrix_indices(:,4),trans_matrix_indices(:,3))];
end

% remove empty rows
state_transition_table=state_transition_table(all(~isnan(state_transition_table),2),:);

% for the kinetic matrix we also need the diagonal elements which are (-1)*sum(rates of outgoing transitions)
nonsink_states=unique(state_transition_table(:,1)); param_table_size=size(transition_rates_table);
% class_val=class(transition_rates_table);
if isa(transition_rates_table,'sym')
    diag_elems_table=sym(zeros(length(nonsink_states),3));
else
    diag_elems_table=zeros(length(nonsink_states),3);
end

n_nonsink_states=length(nonsink_states);

% this loop is the most time consuming
for k=1:length(nonsink_states)
    if rem(k,round(n_nonsink_states/10))==0
        disp(strcat('nonterminal states, ', num2str(round(k/n_nonsink_states,3)) ))
    end
    % which node?
    param_node_ind=state_transition_table(state_transition_table(:,1)==nonsink_states(k),3); 
    % up or down?
    param_updown_ind=state_transition_table(state_transition_table(:,1)==nonsink_states(k),4);
    % diag_elems(nonsink_states(i)) = -sum(transition_rates_table(sub2ind(param_table_size,param_updown_ind,param_node_ind)));
    diag_elems_table(k,:) = [nonsink_states(k), nonsink_states(k), -sum(transition_rates_table(sub2ind(param_table_size,param_updown_ind,param_node_ind)))];
end

disp('done w indexing nonsink states')

% this is index of all nonzero transitions
nondiag_elems_table = [state_transition_table(:,1:2) transition_rates_table(sub2ind(param_table_size,state_transition_table(:,4),state_transition_table(:,3)))];
% this is the full table to generate kinetic matrix K
K_sparse_table = [nondiag_elems_table; diag_elems_table];

disp('done w building K sparse table')

if ~isa(transition_rates_table,'sym')
K_sparse = sparse(K_sparse_table(:,2), K_sparse_table(:,1), K_sparse_table(:,3), size(truth_table_inputs,1), size(truth_table_inputs,1));
A_sparse = transpose(K_sparse/sum(transition_rates_table(:))) + sparse(eye(size(K_sparse)));
else
    % K_sparse=zeros(size(truth_table_inputs,1), size(truth_table_inputs,1)); A_sparse=zeros(size(truth_table_inputs,1), size(truth_table_inputs,1)); 
    K_sparse=sym(zeros(2^numel(nodes),2^numel(nodes))); K_sparse(sub2ind(size(K_sparse),K_sparse_table(:,2),K_sparse_table(:,1)))=K_sparse_table(:,3);
    A_sparse=transpose( K_sparse + eye(size(K_sparse))*sum(transition_rates_table(:)) )/sum(transition_rates_table(:));
   
end

else
   disp('DIMENSION MISMATCH between truth table and node list!!') 
end