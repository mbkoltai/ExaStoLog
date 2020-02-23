function update_matrix =fcn_truthtable(nodes)

n=numel(nodes); list_binary_states=fliplr(rem(floor((0:((2^n)-1)).'*pow2(0:-1:-n+1)),2));

update_matrix = [list_binary_states(:,1),...
list_binary_states(:,2),...
list_binary_states(:,3)|list_binary_states(:,14),...
list_binary_states(:,4)|list_binary_states(:,14),...
list_binary_states(:,14),...
list_binary_states(:,14),...
list_binary_states(:,5)|list_binary_states(:,6),...
(list_binary_states(:,6)|list_binary_states(:,5)|list_binary_states(:,8))&~list_binary_states(:,1),...
list_binary_states(:,8),...
list_binary_states(:,8)&~list_binary_states(:,2),...
list_binary_states(:,8)&(list_binary_states(:,9)|list_binary_states(:,10)),...
list_binary_states(:,12)&~list_binary_states(:,2),...
list_binary_states(:,4)&list_binary_states(:,12),...
~list_binary_states(:,11)&~list_binary_states(:,13)&~list_binary_states(:,3),...
~list_binary_states(:,11)&~list_binary_states(:,3),...
~list_binary_states(:,11)&~list_binary_states(:,13),...
list_binary_states(:,14)&(~list_binary_states(:,11)&~list_binary_states(:,3)),...
~list_binary_states(:,17),...
list_binary_states(:,14),...
~list_binary_states(:,11)&~list_binary_states(:,3),...
(~list_binary_states(:,16)|~list_binary_states(:,15))&~list_binary_states(:,2),...
(list_binary_states(:,18)|list_binary_states(:,21))&~list_binary_states(:,23),...
list_binary_states(:,23)|(list_binary_states(:,19)&list_binary_states(:,20)&~list_binary_states(:,21))];