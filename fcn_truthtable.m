function update_matrix =fcn_truthtable(nodes)

n=numel(nodes); list_binary_states=fliplr(rem(floor((0:((2^n)-1)).'*pow2(0:-1:-n+1)),2));

update_matrix = [~list_binary_states(:,2),...
~list_binary_states(:,1)&list_binary_states(:,3),...
list_binary_states(:,2)|list_binary_states(:,3)];