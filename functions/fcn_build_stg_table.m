function state_transitions_inds=fcn_build_stg_table(truth_table_filename,nodes)

n_nodes=numel(nodes); 

fcn_name_str=strrep(truth_table_filename,'.m',''); fcn_name=str2func(fcn_name_str); 
update_table=sparse(fcn_name(nodes)); 

up_trans_source=arrayfun(@(x) intersect(find(update_table(:,x)),fcn_state_inds(0,n_nodes,x)), 1:size(update_table,2),'un',0);
down_trans_source=arrayfun(@(x) intersect(find(~update_table(:,x)), fcn_state_inds(1,n_nodes,x)), 1:size(update_table,2),'un',0);

down_trans_target=arrayfun(@(x) [[down_trans_source{x}]-2^(x-1) repmat([x 2],numel(down_trans_source{x}),1)], 1:numel(down_trans_source), 'un',0);
up_trans_target=arrayfun(@(x) [[up_trans_source{x}]+2^(x-1) repmat([x 1],numel(up_trans_source{x}),1)], 1:numel(up_trans_source),'un',0);

state_transitions_inds=[ [cell2mat(vertcat(down_trans_source(:))); cell2mat(vertcat(up_trans_source(:)))] ...
           [cell2mat(vertcat(down_trans_target(:))); cell2mat(vertcat(up_trans_target(:)))] ];

% %
% all_binary_states=sparse( rem(floor([0:((2^n_nodes)-1)].'*pow2(0:-1:-n_nodes+1)),2)); % dim_matr=2^n_nodes;
% all_binary_states_decim = sum( bsxfun(@times,all_binary_states,(2.^fliplr(0:n_nodes-1))) ,2);
% 
% % delete(fcn_name)
% if isa(strrep(truth_table_filename,'.m',''),'function_handle')
%     disp('handle exists!!')
%     clear(fcn_name) %  clear(func2str(yourfunctionhandle))
% end
% fcn_name_str=strrep(truth_table_filename,'.m',''); fcn_name=str2func(fcn_name_str); 
% update_table=sparse(fcn_name(nodes)); 
% 
% % if size(update_table,2)~=n_nodes
% %     clear(fcn_name);
% %     update_table=fcn_name(nodes);
% % end
% 
% trans_source_states=cell(n_nodes,1); trans_target_states=trans_source_states; 
% up_down_inds=trans_source_states; node_inds=up_down_inds; % update_asynchr_table=all_binary_states;
% 
% for k=1:n_nodes
%     disp(strcat('building transitions of variable #', num2str(k)));
%     % update_asynchr_table=all_binary_states; update_asynchr_table(:,k)=update_table(:,k);
%     if k==1; x=[]; else; x=all_binary_states(:,1:max([1 k-1])); end
%     if k==n_nodes; y=[]; else; y=all_binary_states(:,min(k+1,n_nodes):n_nodes); end
%     update_asynchr_table=[x update_table(:,k) y];
%     trans_source_states{k}=all_binary_states(:,k)~=update_asynchr_table(:,k); 
%     % this is the time consuming step
%     new_states = update_asynchr_table(trans_source_states{k},:);
%     [~,I_outputs,I_inputs]=intersect(all_binary_states_decim, ... % (~trans_source_states{k})
%                                       sum( bsxfun(@times,new_states,(2.^fliplr(0:n_nodes-1)) ),2) );
%     % 
%     [~,b]=ismember(1:numel(I_inputs),I_inputs); % [~,b]=intersect(I_inputs,(1:numel(I_inputs))'); 
%     trans_target_states{k} = I_outputs(b);
%     up_down_inds{k} = full(sum(new_states - all_binary_states(trans_source_states{k},:),2 ));
%     node_inds{k} = repelem(k,numel(I_inputs),1);
% end
% 
% % convert -1 to 2 (down rates)
% up_down_inds_mat=cell2mat(up_down_inds); trans_source_states_mat=cell2mat(cellfun(@(x) find(x),trans_source_states,'un',0));
% trans_target_states_mat=cell2mat(trans_target_states);
% up_down_inds_arr = sign(up_down_inds_mat).^2 + (up_down_inds_mat<0);
% 
% % [source_state target_state node_ind up_down_ind]
% state_transitions_inds=[trans_source_states_mat trans_target_states_mat cell2mat(node_inds) up_down_inds_arr];
