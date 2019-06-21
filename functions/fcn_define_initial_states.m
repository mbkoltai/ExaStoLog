function x0=fcn_define_initial_states(initial_on_nodes,dom_prob,nodes)

n_nodes=numel(nodes); truth_table_inputs=rem(floor([0:((2^n_nodes)-1)].'*pow2(0:-1:-n_nodes+1)),2);
% define initial values
x0=zeros(1,2^n_nodes)'; 
% defining a dominant initial state (eg. dom_prob=0.8, ie. 80% probability)
initial_on_nodes_inds=find(sum(cell2mat(arrayfun(@(x) strcmp(nodes,initial_on_nodes{x} ), 1:numel(initial_on_nodes),'un',0)')));
% dom_prob=0.8; 
initial_state=zeros(1,numel(nodes)); initial_state(initial_on_nodes_inds)=1;
x0(ismember(truth_table_inputs,initial_state,'rows'))=dom_prob; 

% take those states where selected variables have a value of 1, and we want them to have a nonzero probability
sel_states=all(truth_table_inputs(:,initial_on_nodes_inds)');
% create a vector of random probabilities for these states, with a sum of (1-dom_prob)
rand_vect=abs(rand(sum(sel_states)-1,1)); rand_vect=rand_vect/( sum(rand_vect)/(1-dom_prob) );
x0(~ismember(truth_table_inputs,initial_state,'rows') & sel_states') = rand_vect;
