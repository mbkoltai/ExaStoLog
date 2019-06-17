function B_state_transitions_inds=fcn_stg_table_subgraph(stg_table,cell_subgraphs,counter)

B_state_transitions_inds=stg_table(ismember(stg_table(:,1),cell_subgraphs{counter}),:);
[~,b]=ismember(B_state_transitions_inds(:,1),cell_subgraphs{counter}); [~,d]=ismember(B_state_transitions_inds(:,2),cell_subgraphs{counter}); 
B_state_transitions_inds(:,1)=b; B_state_transitions_inds(:,2)=d; 
